=head1 LICENSE

Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy self the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, sselftware
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

# EnsEMBL module for Bio::EnsEMBL::VEP::OutputFactory
#
#

=head1 NAME

Bio::EnsEMBL::VEP::OutputFactory - base class for generating VEP output

=head1 SYNOPSIS

# actually gets a Bio::EnsEMBL::VEP::OutputFactory::VCF
my $of = Bio::EnsEMBL::VEP::OutputFactory->new({
  config => $config,
  format => 'vcf'
});

# print headers
print "$_\n" for @{$of->headers};

# print output
print "$_\n" for @{$of->get_all_lines_by_InputBuffer($ib)};

=head1 DESCRIPTION

The OutputFactory class is a base class used to generate VEP output.

It should not be invoked directly, but contains the bulk of the methods
used to transform the Ensembl API objects into "flat" structures suitable
for writing to output.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::OutputFactory;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Scalar::Util qw(looks_like_number);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Constants;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::VEP::Utils qw(format_coords merge_arrays);
use Bio::EnsEMBL::VEP::Constants;

use Bio::EnsEMBL::VEP::OutputFactory::VEP_output;
use Bio::EnsEMBL::VEP::OutputFactory::VCF;
use Bio::EnsEMBL::VEP::OutputFactory::Tab;

our $CAN_USE_JSON;

BEGIN {
  if(eval q{ use Bio::EnsEMBL::VEP::OutputFactory::JSON; 1 }) {
    $CAN_USE_JSON = 1;
  }
}

my %SO_RANKS = map {$_->SO_term => $_->rank} values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;

my %FORMAT_MAP = (
  'vcf'     => 'VCF',
  'ensembl' => 'VEP_output',
  'vep'     => 'VEP_output',
  'tab'     => 'Tab',
  'json'    => 'JSON',
);

my %DISTANCE_CONS = (upstream_gene_variant => 1, downstream_gene_variant => 1);

my %FREQUENCY_KEYS = (
  af        => ['AF'],
  af_1kg    => [qw(AFR AMR ASN EAS EUR SAS)],
  af_esp    => [qw(AA EA)],
  af_exac   => [('ExAC', map {'ExAC_'.$_} qw(Adj AFR AMR EAS FIN NFE OTH SAS))],
  af_gnomad => [('gnomAD', map {'gnomAD_'.$_} qw(AFR AMR ASJ EAS FIN NFE OTH SAS))],
);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config,
                 format => string (vcf, ensembl, vep, tab, json)
               }
  Example    : $of = Bio::EnsEMBL::VEP::OutputFactory({
                 config => $config,
                 format => 'tab'
               });
  Description: Create a new Bio::EnsEMBL::VEP::OutputFactory object.
               Will return the sub-class determined by the format arg.
  Returntype : Bio::EnsEMBL::VEP::OutputFactory
  Exceptions : throws on invalid format or if JSON format
               requested and JSON module not installed
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    no_stats
    no_intergenic
    process_ref_homs
    coding_only
    terms
    output_format
    no_escape
    pick_order
    allele_number
    show_ref_allele
    use_transcript_ref
    vcf_string

    most_severe
    summary
    per_gene
    pick
    pick_allele
    pick_allele_gene
    flag_pick
    flag_pick_allele
    flag_pick_allele_gene
    
    variant_class
    af
    af_1kg
    af_esp
    af_exac
    af_gnomad
    max_af
    pubmed
    clin_sig_allele

    numbers
    domains
    symbol
    ccds
    xref_refseq
    refseq
    merged
    protein
    uniprot
    canonical
    biotype
    mane
    mane_select
    mane_plus_clinical
    tsl
    appris
    transcript_version
    gene_phenotype
    mirna
    ambiguity
    var_synonyms

    total_length
    hgvsc
    hgvsp
    hgvsg
    hgvsg_use_accession
    spdi
    sift
    polyphen
    polyphen_analysis

    cell_type
    shift_3prime
    shift_genomic
    remove_hgvsp_version
  )]);

  my $hashref = $_[0];

  if(my $format = $hashref->{format}) {

    delete $hashref->{format};

    $format = lc($format);
    throw("ERROR: Unknown or unsupported output format $format\n") unless $FORMAT_MAP{$format};


    if($FORMAT_MAP{$format} eq 'JSON') {
      throw("ERROR: Cannot use format $format without JSON module installed\n") unless $CAN_USE_JSON;
    }

    my $class = 'Bio::EnsEMBL::VEP::OutputFactory::'.$FORMAT_MAP{$format};
    return $class->new({%$hashref, config => $self->config});
  }

  $self->{header_info} = $hashref->{header_info} if $hashref->{header_info};

  $self->{plugins} = $hashref->{plugins} if $hashref->{plugins};;

  return $self;
}


### Top level methods
#####################


=head2 get_all_lines_by_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Example    : $lines = $of->get_all_lines_by_InputBuffer($ib);
  Description: Gets all lines (strings suitable for writing to output) given
               an annotated input buffer.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub get_all_lines_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  # output_hash_to_line will be implemented in the child class as its behaviour varies by output type
  return [map {$self->output_hash_to_line($_)} @{$self->get_all_output_hashes_by_InputBuffer($buffer)}];
}


=head2 get_all_output_hashes_by_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Example    : $hashes = $of->get_all_output_hashes_by_InputBuffer($ib);
  Description: Gets all hashrefs of data (suitable for serialised writing e.g. JSON) given
               an annotated input buffer.
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub get_all_output_hashes_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  # rejoin variants that have been split
  # this can happen when using --minimal
  $self->rejoin_variants_in_InputBuffer($buffer) if $buffer->rejoin_required;


  map {@{$self->reset_shifted_positions($_)}}
    @{$buffer->buffer};
  
  return [
    map {@{$self->get_all_output_hashes_by_VariationFeature($_)}}
    @{$buffer->buffer}
  ];
}

sub reset_shifted_positions {
  my $self = shift;
  my $vf = shift;
  return [] if ref($vf) eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature';
  my @tvs = $vf->get_all_TranscriptVariations();
  foreach my $tv (@{$tvs[0]})
  {
      map { bless $_, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele' } 
          @{ $tv->get_all_BaseVariationFeatureOverlapAlleles };
        
      ## Obtains all relevant $tva objects and removes shifted positions from 
      ## CDS CDNA and Protein positions   
      map { $_->clear_shifting_variables} 
          @{ $tv->get_all_BaseVariationFeatureOverlapAlleles };
  }

  return \@tvs;
}



=head2 get_all_output_hashes_by_VariationFeature

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : $hashes = $of->get_all_output_hashes_by_InputBuffer($vf);
  Description: This wrapper method does the following:
                1) fetches all the VariationFeatureOverlapAllele objects
                2) filters them if any filtering options are applied
                3) converts them into a flat hashref with the relevant output data
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : get_all_output_hashes_by_InputBuffer()
  Status     : Stable

=cut

sub get_all_output_hashes_by_VariationFeature {
  my $self = shift;
  my $vf = shift;

  # get the basic initial hash; this basically contains the location and ID
  my $hash = $self->VariationFeature_to_output_hash($vf);

  return $self->get_all_VariationFeatureOverlapAllele_output_hashes($vf, $hash);
}


=head2 get_all_VariationFeatureOverlapAllele_output_hashes

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg 2      : hashref of data obtained from VariationFeature
  Example    : $hashes = $of->get_all_VariationFeatureOverlapAllele_output_hashes($vf, $hash);
  Description: Takes the initial hashref of data obtained from just
               the VariationFeature (i.e. location-based data only)
               and then creates a hashref expanding on this for each
               VariationFeatureOverlapAllele (e.g. TranscriptVariationAllele)
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : get_all_output_hashes_by_VariationFeature()
  Status     : Stable

=cut

sub get_all_VariationFeatureOverlapAllele_output_hashes {
  my $self = shift;
  my $vf = shift;
  my $hash = shift;

  my @return;

  # we want to handle StructuralVariationFeatures differently
  my $method =  sprintf(
    'get_all_%sOverlapAlleles',
    ref($vf) eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature'
      ? 'StructuralVariation'
      : 'VariationFeature'
  );

  my $vfoas = $self->$method($vf);

  # summary, most_severe don't need most of the downstream logic from hereon
  return $self->summary_only($vf, $hash, $vfoas) if $self->{summary} || $self->{most_severe};

  foreach my $vfoa(@$vfoas) {

    # copy the initial VF-based hash so we're not overwriting
    my %copy = %$hash;

    # we have a method defined for each sub-class of VariationFeatureOverlapAllele
    my $method = (split('::', ref($vfoa)))[-1].'_to_output_hash';
    my $output = $self->$method($vfoa, \%copy, $vf);

    # run plugins
    $output = $self->run_plugins($vfoa, $output, $vf);

    # log stats
    $self->stats->log_VariationFeatureOverlapAllele($vfoa, $output) unless $self->{no_stats};

    push @return, $output if $output;
  }

  return \@return;
}


=head2 summary_only

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg 2      : hashref of data obtained from VariationFeature
  Arg 3      : arrayref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Example    : $hashes = $of->summary_only($vf, $hash, $vfoas);
  Description: Squashes data from a listref of VariationFeatureOverlapAlleles
               into one hashref, using either a summary (list all consequence types)
               or most severe method.
  Returntype : arrayref of one hashref (for compliance with other methods)
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub summary_only {
  my ($self, $vf, $hash, $vfoas) = @_;

  my $term_method = $self->{terms}.'_term';

  my @ocs = sort {$a->rank <=> $b->rank} map {@{$_->get_all_OverlapConsequences}} @$vfoas;

  if(@ocs) {

    # summary is just all unique consequence terms
    if($self->{summary}) {
      my (@cons, %seen);
      foreach my $con(map {$_->$term_method || $_->SO_term} @ocs) {
        push @cons, $con unless $seen{$con};
        $seen{$con} = 1;
      }
      $hash->{Consequence} = \@cons;
    }

    # most severe is the consequence term with the lowest rank
    else {
      $hash->{Consequence} = [$ocs[0]->$term_method || $ocs[0]->SO_term];
    }

    # unless(defined($config->{no_stats})) {
    #   $config->{stats}->{consequences}->{$_}++ for split(',', $hash->{Consequence});
    # }
  }
  else {
    $self->warning_msg("Unable to assign consequence type");
  }

  return [$hash];
}


=head2 get_all_VariationFeatureOverlapAlleles

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : $vfoas = $of->get_all_VariationFeatureOverlapAlleles($vf);
  Description: Gets all VariationFeatureOverlapAllele objects for given
               VariationFeature, applying any filtering as configured.
               This might be restricting to coding only, or filtering using
               one of the pick* functions.
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub get_all_VariationFeatureOverlapAlleles {
  my $self = shift;
  my $vf = shift;

  # no intergenic?
  return [] if $self->{no_intergenic} && defined($vf->{intergenic_variation});

  # get all VFOAs
  # need to be sensitive to whether --coding_only is switched on
  my $vfos = $vf->get_all_VariationFeatureOverlaps;

  # if coding only filter TranscriptVariations down to coding ones
  # leave others intact, otherwise doesn't make sense to do --regulatory and --coding_only
  if($self->{coding_only}) {
    my @new;

    foreach my $vfo(@$vfos) {
      if($vfo->isa('Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap')) {
        push @new, $vfo if $vfo->can('affects_cds') && $vfo->affects_cds;
      }
      else {
        push @new, $vfo;
      }
    }

    $vfos = \@new;
  }

  # method name stub for getting *VariationAlleles
  my $allele_method = $self->{process_ref_homs} ? 'get_all_' : 'get_all_alternate_';  
  my $method = $allele_method.'VariationFeatureOverlapAlleles';

  return $self->filter_VariationFeatureOverlapAlleles([map {@{$_->$method}} @{$vfos}]);
}


=head2 get_all_StructuralVariationOverlapAlleles

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : $svos = $of->get_all_StructuralVariationOverlapAlleles($vf);
  Description: Gets all StructuralVariationOverlap objects for given
               StructuralVariation, applying any filtering as configured.
               This might be restricting to coding only, or filtering using
               one of the pick* functions.
  Returntype : arrayref of Bio::EnsEMBL::Variation::StructuralVariationOverlap
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub get_all_StructuralVariationOverlapAlleles {
  my $self = shift;
  my $svf = shift;

  # no intergenic?
  my $isv = $svf->get_IntergenicStructuralVariation(1);
  return [] if $self->{no_intergenic} && $isv;

  # get all VFOAs
  # need to be sensitive to whether --coding_only is switched on
  my $vfos;

  # if coding only just get transcript & intergenic ones
  if($self->{coding_only}) {
    @$vfos = grep {defined($_)} (
      @{$svf->get_all_TranscriptStructuralVariations},
      $isv
    );
  }
  else {
    $vfos = $svf->get_all_StructuralVariationOverlaps;
  }

  # grep out non-coding?
  @$vfos = grep {$_->can('affects_cds') && $_->affects_cds} @$vfos if $self->{coding_only};

  return $self->filter_StructuralVariationOverlapAlleles([map {@{$_->get_all_StructuralVariationOverlapAlleles}} @{$vfos}]);
}


=head2 filter_VariationFeatureOverlapAlleles

  Arg 1      : arrayref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Example    : $filtered = $of->filter_VariationFeatureOverlapAlleles($vfoas);
  Description: Filters VariationFeatureOverlapAllele objects by various criteria,
               choosing one per variant, allele, gene. It can also flag instead
               of filtering.
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAlleles()
  Status     : Stable

=cut

sub filter_VariationFeatureOverlapAlleles {
  my $self = shift;
  my $vfoas = shift;

  return [] unless $vfoas && scalar @$vfoas;

  # pick worst?
  if($self->{pick}) {
    return [$self->pick_worst_VariationFeatureOverlapAllele($vfoas)];
  }

  # pick worst per allele?
  elsif($self->{pick_allele}) {
    my %by_allele;
    push @{$by_allele{$_->variation_feature_seq}}, $_ for @$vfoas;
    return [map {$self->pick_worst_VariationFeatureOverlapAllele($by_allele{$_})} keys %by_allele];
  }

  # pick per gene?
  elsif($self->{per_gene}) {
    return $self->pick_VariationFeatureOverlapAllele_per_gene($vfoas);
  }

  # pick worst per allele and gene?
  elsif($self->{pick_allele_gene}) {
    my %by_allele;
    push @{$by_allele{$_->variation_feature_seq}}, $_ for @$vfoas;
    return [map {@{$self->pick_VariationFeatureOverlapAllele_per_gene($by_allele{$_})}} keys %by_allele];
  }

  # flag picked?
  elsif($self->{flag_pick}) {
    if(my $worst = $self->pick_worst_VariationFeatureOverlapAllele($vfoas)) {
      $worst->{PICK} = 1;
    }
  }

  # flag worst per allele?
  elsif($self->{flag_pick_allele}) {
    my %by_allele;
    push @{$by_allele{$_->variation_feature_seq}}, $_ for @$vfoas;
    $self->pick_worst_VariationFeatureOverlapAllele($by_allele{$_})->{PICK} = 1 for keys %by_allele;
  }

  # flag worst per allele and gene?
  elsif($self->{flag_pick_allele_gene}) {
    my %by_allele;
    push @{$by_allele{$_->variation_feature_seq}}, $_ for @$vfoas;
    map {$_->{PICK} = 1} map {@{$self->pick_VariationFeatureOverlapAllele_per_gene($by_allele{$_})}} keys %by_allele;
  }

  return $vfoas;
}


=head2 filter_StructuralVariationOverlapAlleles

  Arg 1      : arrayref of Bio::EnsEMBL::Variation::StructuralVariationOverlapAlleles
  Example    : $filtered = $of->filter_StructuralVariationOverlapAlleles($svoas);
  Description: Filters StructuralVariationOverlapAlleles objects by various criteria,
               choosing one per variant, allele, gene. It can also flag instead
               of filtering.
  Returntype : arrayref of Bio::EnsEMBL::Variation::StructuralVariationOverlapAlleles
  Exceptions : none
  Caller     : get_all_StructuralVariationOverlapAlleles()
  Status     : Stable

=cut

sub filter_StructuralVariationOverlapAlleles {
  my $self = shift;
  my $svoas = shift;

  return [] unless $svoas && scalar @$svoas;

  # pick worst? pick worst per allele?
  if($self->{pick} || $self->{pick_allele}) {
    return [$self->pick_worst_VariationFeatureOverlapAllele($svoas)];
  }

  # pick per gene? pick worst per allele and gene?
  elsif($self->{per_gene} || $self->{pick_allele_gene}) {
    return $self->pick_VariationFeatureOverlapAllele_per_gene($svoas);
  }

  # flag picked? flag worst per allele?
  elsif($self->{flag_pick} || $self->{flag_pick_allele}) {
    if(my $worst = $self->pick_worst_VariationFeatureOverlapAllele($svoas)) {
      $worst->{PICK} = 1;
    }
  }

  # flag worst per allele and gene?
  elsif($self->{flag_pick_allele_gene}) {
    map {$_->{PICK} = 1} @{$self->pick_VariationFeatureOverlapAllele_per_gene($svoas)};
  }

  return $svoas;

}


=head2 pick_worst_VariationFeatureOverlapAllele

  Arg 1      : arrayref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Example    : $picked = $of->pick_worst_VariationFeatureOverlapAllele($vfoas);
  Description: Selects one VariationFeatureOverlapAllele from a list using criteria
               defined in the param pick_order. Criteria are in this default order:
                1: canonical
                2: transcript support level
                3: biotype (protein coding favoured)
                4: consequence rank
                5: transcript length
                6: transcript from Ensembl?
                7: transcript from RefSeq?
  Returntype : Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Exceptions : none
  Caller     : filter_VariationFeatureOverlapAlleles(),
               filter_StructuralVariationOverlapAlleles()
  Status     : Stable

=cut

sub pick_worst_VariationFeatureOverlapAllele {
  my $self = shift;
  my $vfoas = shift || [];

  my @vfoa_info;

  return $vfoas->[0] if scalar @$vfoas == 1;

  foreach my $vfoa(@$vfoas) {

    # create a hash self info for this VFOA that will be used to rank it
    my $info = {
      vfoa => $vfoa,
      rank => undef,

      # these will only be used by transcript types, default to 1 for others
      # to avoid writing an else clause below
      mane => 1,
      canonical => 1,
      ccds => 1,
      length => 0,
      biotype => 1,
      tsl => 100,
      appris => 100,
      ensembl => 1,
      refseq => 1,
    };

    if($vfoa->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele')) {
      my $tr = $vfoa->feature;

      # 0 is "best"
      $info->{mane} = scalar(grep {$_->code eq 'MANE_Select'}  @{$tr->get_all_Attributes()}) ? 0 : 1;
      $info->{canonical} = $tr->is_canonical ? 0 : 1;
      $info->{biotype} = $tr->biotype eq 'protein_coding' ? 0 : 1;
      $info->{ccds} = $tr->{_ccds} && $tr->{_ccds} ne '-' ? 0 : 1;
      $info->{lc($tr->{_source_cache})} = 0 if exists($tr->{_source_cache});

      # "invert" length so longer is best
      $info->{length} = 0 - (
        $tr->translation ?
        length($tr->{_variation_effect_feature_cache}->{translateable_seq} || $tr->translateable_seq) :
        $tr->length()
      );

      # lower TSL is best
      if(my ($tsl) = @{$tr->get_all_Attributes('TSL')}) {
        if($tsl->value =~ m/tsl(\d+)/) {
          $info->{tsl} = $1 if $1;
        }
      }

      # lower APPRIS is best
      if(my ($appris) = @{$tr->get_all_Attributes('appris')}) {
        if($appris->value =~ m/([A-Za-z]).+(\d+)/) {
          my ($type, $grade) = ($1, $2);

          # values are principal1, principal2, ..., alternative1, alternative2
          # so add 10 to grade if alternative
          $grade += 10 if substr($type, 0, 1) eq 'a';

          $info->{appris} = $grade if $grade;
        }
      }
    }

    push @vfoa_info, $info;
  }
  if(scalar @vfoa_info) {
    my @order = @{$self->{pick_order}};
    my $picked;

    # go through each category in order
    foreach my $cat(@order) {

      # get ranks here as it saves time
      if($cat eq 'rank') {
        foreach my $info(@vfoa_info) {
          my @ocs = sort {$a->rank <=> $b->rank} @{$info->{vfoa}->get_all_OverlapConsequences};
          $info->{rank} = scalar @ocs ? $SO_RANKS{$ocs[0]->SO_term} : 1000;
        }
      }

      # sort on that category
      @vfoa_info = sort {$a->{$cat} <=> $b->{$cat}} @vfoa_info;

      # take the first (will have the lowest value self $cat)
      $picked = shift @vfoa_info;
      my @tmp = ($picked);

      # now add to @tmp those vfoas that have the same value self $cat as $picked
      push @tmp, shift @vfoa_info while @vfoa_info && $vfoa_info[0]->{$cat} eq $picked->{$cat};

      # if there was only one, return
      return $picked->{vfoa} if scalar @tmp == 1;

      # otherwise shrink the array to just those that had the lowest
      # this gives fewer to sort on the next round
      @vfoa_info = @tmp;
    }

    # probably shouldn't get here, but if we do, return the first
    return $vfoa_info[0]->{vfoa};
  }

  return undef;
}


=head2 pick_VariationFeatureOverlapAllele_per_gene

  Arg 1      : arrayref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Example    : $picked = $of->pick_VariationFeatureOverlapAllele_per_gene($vfoas);
  Description: Selects one VariationFeatureOverlapAllele per gene affected by the
               given list, using pick_worst_VariationFeatureOverlapAllele to select
               within in each gene.
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele
  Exceptions : none
  Caller     : filter_VariationFeatureOverlapAlleles(),
               filter_StructuralVariationOverlapAlleles()
  Status     : Stable

=cut

sub pick_VariationFeatureOverlapAllele_per_gene {
  my $self = shift;
  my $vfoas = shift;

  my @return;
  my @tvas;

  # pick out TVAs
  foreach my $vfoa(@$vfoas) {
    if($vfoa->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele')) {
      push @tvas, $vfoa;
    }
    else {
      push @return, $vfoa;
    }
  }

  # sort the TVA objects into a hash by gene
  my %by_gene;

  foreach my $tva(@tvas) {
    my $gene = $tva->transcript->{_gene_stable_id} || $self->{ga}->fetch_by_transcript_stable_id($tva->transcript->stable_id)->stable_id;
    push @{$by_gene{$gene}}, $tva;
  }

  foreach my $gene(keys %by_gene) {
    push @return, grep {defined($_)} $self->pick_worst_VariationFeatureOverlapAllele($by_gene{$gene});
  }

  return \@return;
}


### Transforming methods
### Typically these methods will transform an API object into hash-able data
### Most will take the object to be transformed and the hash into which the
### data are inserted as the arguments, and return the hashref
############################################################################


=head2 VariationFeature_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature
  Example    : $hashref = $of->VariationFeature_to_output_hash($vf);
  Description: "Flattens" a VariationFeature to a simple hash, adding information
               pertaining just to the variant's genomic location.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_output_hashes_by_VariationFeature()
  Status     : Stable

=cut

sub VariationFeature_to_output_hash {
  my $self = shift;
  my $vf = shift;

  my $hash = {
    Uploaded_variation  => $vf->variation_name ne '.' ? $vf->variation_name : ($vf->{original_chr} || $vf->{chr}).'_'.$vf->{start}.'_'.($vf->{allele_string} || $vf->{class_SO_term}),
    Location            => ($vf->{chr} || $vf->seq_region_name).':'.format_coords($vf->{start}, $vf->{end}),
  };

  my $converted_to_vcf = $vf->to_VCF_record;

  my $alt_allele_vcf = ${$converted_to_vcf}[4];

  if($self->{vcf_string} || (defined($self->{_config}->{_params}->{fields}) && grep(/vcf_string/, @{$self->{_config}->{_params}->{fields}}))){
    if($alt_allele_vcf =~ /,/){
      my @list_vcfs;
      my @alt_splited_list = split(q(,), $alt_allele_vcf);

      foreach my $alt_splited (@alt_splited_list){
        push(@list_vcfs, $vf->{chr}.'-'.${$converted_to_vcf}[1].'-'.${$converted_to_vcf}[3].'-'.$alt_splited);
      }
      $hash->{vcf_string} = \@list_vcfs;
    }
    else{
      $hash->{vcf_string} = $vf->{chr}.'-'.${$converted_to_vcf}[1].'-'.${$converted_to_vcf}[3].'-'.${$converted_to_vcf}[4];
    }
  }

  # get variation synonyms for Variant Recoder
  # if the tool is Variant Recoder fetches the variation synonyms from the database
  if($self->{_config}->{_params}->{is_vr} && $self->{var_synonyms}){
    my $variation = $vf->variation();
    my $var_synonyms = $variation->get_all_synonyms('', 1);
    $hash->{var_synonyms} = $var_synonyms;
  }

  # overlapping SVs
  if($vf->{overlapping_svs}) {
    $hash->{SV} = [sort keys %{$vf->{overlapping_svs}}];
  }

  # variant class
  $hash->{VARIANT_CLASS} = $vf->class_SO_term() if $self->{variant_class};

  # individual?
  if(defined($vf->{individual})) {
    $hash->{IND} = $vf->{individual};

    # zygosity
    if(defined($vf->{genotype})) {
    my %unique = map {$_ => 1} @{$vf->{genotype}};
    $hash->{ZYG} = (scalar keys %unique > 1 ? 'HET' : 'HOM').(defined($vf->{hom_ref}) ? 'REF' : '');
    }
  }

  # minimised?
  $hash->{MINIMISED} = 1 if $vf->{minimised};
  
  if(ref($vf) eq 'Bio::EnsEMBL::Variation::VariationFeature') {
    my $ambiguity_code = $vf->ambig_code();
    
    if($self->{ambiguity} && defined($ambiguity_code)) {
      $hash->{AMBIGUITY} = $ambiguity_code;
    }
  }

  # custom annotations
  foreach my $custom_name(keys %{$vf->{_custom_annotations} || {}}) {
    $self->_add_custom_annotations_to_hash(
      $hash,
      $custom_name,
      [
        grep {!exists($_->{allele})}
        @{$vf->{_custom_annotations}->{$custom_name}}
      ]
    );
  }

  # nearest
  $hash->{NEAREST} = $vf->{nearest} if $vf->{nearest};

  # check_ref tests
  $hash->{CHECK_REF} = 'failed' if defined($vf->{check_ref_failed});

  $self->stats->log_VariationFeature($vf, $hash) unless $self->{no_stats};

  return $hash;
}


=head2 add_colocated_variant_info

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->add_colocated_variant_info($vf, $vf_hash);
  Description: Adds co-located variant information to hash
  Returntype : hashref
  Exceptions : none
  Caller     : VariationFeature_to_output_hash()
  Status     : Stable

=cut

sub add_colocated_variant_info {
  my $self = shift;
  my $vf = shift;
  my $hash = shift;

  return unless $vf->{existing} && scalar @{$vf->{existing}};

  my $this_allele = $hash->{Allele};

  my $shifted_allele = $vf->{shifted_allele_string};
  $shifted_allele ||= "";
  my $tmp = {};

  my $clin_sig_allele_exists = 0;
  # use these to sort variants
  my %prefix_ranks = (
    'rs' => 1, # dbSNP

    'cm' => 2, # HGMD
    'ci' => 2,
    'cd' => 2,

    'co' => 3, # COSMIC
  );

  my %clin_sigs;
  
  foreach my $ex(
    sort {
      ($a->{somatic} || 0) <=> ($b->{somatic} || 0) ||

      ($prefix_ranks{lc(substr($a->{variation_name}, 0, 2))} || 100)
      <=>
      ($prefix_ranks{lc(substr($b->{variation_name}, 0, 2))} || 100)
    }
    @{$vf->{existing}}
  ) {

    # check allele match
    if(my $matched = $ex->{matched_alleles}) {
      next unless (grep {$_->{a_allele} eq $this_allele} @$matched) || (grep {$_->{a_allele} eq $shifted_allele} @$matched) ;
    }

    # ID
    push @{$hash->{Existing_variation}}, $ex->{variation_name} if $ex->{variation_name};

    # Variation Synonyms
    # VEP fetches the variation synonyms from the cache
    push @{$hash->{VAR_SYNONYMS}}, $ex->{var_synonyms} if $self->{var_synonyms} && $ex->{var_synonyms} && !$self->{_config}->{_params}->{is_vr};

    # Find allele specific clin_sig data if it exists
    if(defined($ex->{clin_sig_allele}) && $self->{clin_sig_allele} )
    {
      my %cs_hash;
      my @clin_sig_array = split(';', $ex->{clin_sig_allele});       
      foreach my $cs(@clin_sig_array){
        my @cs_split = split(':', $cs);
        $cs_hash{$cs_split[0]} = '' if !defined($cs_hash{$cs_split[0]});
        $cs_hash{$cs_split[0]} .= ',' if $cs_hash{$cs_split[0]} ne ''; 
        $cs_hash{$cs_split[0]} .= $cs_split[1];
      }

      my $hash_ref = \%cs_hash;
      $clin_sigs{$hash_ref->{$this_allele}} = 1 if defined($hash_ref->{$this_allele});
      $clin_sig_allele_exists = 1;
    }

    # clin sig and pubmed?
    push @{$tmp->{CLIN_SIG}}, split(',', $ex->{clin_sig}) if $ex->{clin_sig} && !$clin_sig_allele_exists;
    push @{$tmp->{PUBMED}}, split(',', $ex->{pubmed}) if $self->{pubmed} && $ex->{pubmed};

    # somatic?
    push @{$tmp->{SOMATIC}}, $ex->{somatic} ? 1 : 0;

    # phenotype or disease
    push @{$tmp->{PHENO}}, $ex->{phenotype_or_disease} ? 1 : 0;   
  }

  # post-process to remove all-0, e.g. SOMATIC
  foreach my $key(keys %$tmp) {
    delete $tmp->{$key} unless grep {$_} @{$tmp->{$key}};
  }
  
  # post-process to merge var synonyms into one entry so we can control the delimiter
  $hash->{VAR_SYNONYMS} = join '--', @{$hash->{VAR_SYNONYMS}} if defined($hash->{VAR_SYNONYMS});

  my @keys = keys(%clin_sigs);
  $tmp->{CLIN_SIG} = join(';', @keys) if scalar(@keys) && $self->{clin_sig_allele};
 
  # copy to hash
  $hash->{$_} = $tmp->{$_} for keys %$tmp;
  # frequencies used to filter will appear here
  if($vf->{_freq_check_freqs}) {
    my @freqs;

    foreach my $p(keys %{$vf->{_freq_check_freqs}}) {
      foreach my $a(keys %{$vf->{_freq_check_freqs}->{$p}}) {
        push @freqs, sprintf(
          '%s:%s:%s',
          $p,
          $a,
          $vf->{_freq_check_freqs}->{$p}->{$a}
        )
      }
    }

    $hash->{FREQS} = \@freqs;
  }

  return $hash;
}



=head2 add_colocated_frequency_data

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg 2      : hashref $vf_hash
  Arg 3      : hashref $existing_variant_hash
  Example    : $hashref = $of->add_colocated_frequency_data($vf, $vf_hash, $ex);
  Description: Adds co-located variant frequency data information to hash
  Returntype : hashref
  Exceptions : none
  Caller     : VariationFeatureOverlapAllele_to_output_hash()
  Status     : Stable

=cut

sub add_colocated_frequency_data {
  my $self = shift;
  my ($vf, $hash, $ex, $shift_hash) = @_;

  return $hash unless grep {$self->{$_}} keys %FREQUENCY_KEYS or $self->{max_af};

  my @ex_alleles = split('/', $ex->{allele_string});

  # gmaf stored a bit differently, but we can get it in the same format
  $ex->{AF} = $ex->{minor_allele}.':'.$ex->{minor_allele_freq} if $ex->{minor_allele};

  my @keys = keys %FREQUENCY_KEYS;
  @keys = grep {$self->{$_}} @keys unless $self->{max_af};
  
  my $this_allele = $hash->{Allele} ||= '-'; #if exists($hash->{Allele});
  my $this_allele_unshifted = $shift_hash->{alt_orig_allele_string} if defined($shift_hash);
  $this_allele_unshifted ||= "";
  
  my ($matched_allele) = grep {$_->{a_allele} eq $this_allele || $_->{a_allele} eq $this_allele_unshifted} @{$ex->{matched_alleles} || []};

  return $hash unless $matched_allele || (grep {$_ eq 'af'} @keys);

  my $max_af = 0;
  my @max_af_pops;

  foreach my $group(@keys) {
    foreach my $key(grep {$ex->{$_}} @{$FREQUENCY_KEYS{$group}}) {

      my %freq_data;
      my $total = 0;

      # use this to log which alleles we've explicitly seen freqs for
      my %remaining = map {$_ => 1} @ex_alleles;

      # get the frequencies for each allele into a hashref
      foreach my $pair(split(',', $ex->{$key})) {
        my ($a, $f) = split(':', $pair);
        $freq_data{$a} = $f;
        $total += $f;
        delete $remaining{$a} if $remaining{$a};
      }

      # interpolate the frequency for the remaining allele if there's only 1
      # we can only do this reliably for the AF key as only the minor AF is stored
      # for others we expect all ALTs to have a store frequency, those without cannot be reliably interpolated
      my $interpolated = 0;
      if(scalar @ex_alleles == 2 && scalar keys %remaining == 1 && $key eq 'AF') {
        $freq_data{(keys %remaining)[0]} = 1 - $total;
        $interpolated = 1;
      }

      if(
        ($matched_allele && exists($freq_data{$matched_allele->{b_allele}})) ||
        ($interpolated && $freq_data{$this_allele})
      ) {

        my $f =
          $matched_allele && exists($freq_data{$matched_allele->{b_allele}}) ?
          $freq_data{$matched_allele->{b_allele}} :
          $freq_data{$this_allele};

        # record the frequency if requested
        if($self->{$group}) {
          my $out_key = $key eq 'AF' ? 'AF' : $key.'_AF';
          merge_arrays($hash->{$out_key} ||= [], [$f]);
        }

        # update max_af data if required
        # make sure we don't include any combined-level pops
        if($self->{max_af} && $key ne 'AF' && $key ne 'ExAC' && $key ne 'ExAC_Adj' && $key ne 'gnomAD') {
          if($f > $max_af) {
            $max_af = $f;
            @max_af_pops = ($key);
          }
          elsif($f == $max_af) {
            push @max_af_pops, $key;
          }
        }
      }
    }
  }

  # add/update max_af info
  if($self->{max_af} && @max_af_pops) {
    my $current_max = $hash->{MAX_AF} ||= 0;

    if($max_af > $current_max) {
      $hash->{MAX_AF} = $max_af;
      $hash->{MAX_AF_POPS} = [];
    }
    
    push @{$hash->{MAX_AF_POPS}}, @max_af_pops if $max_af >= $current_max;
  }

  delete $ex->{AF};

  return $hash;
}


=head2 _add_custom_annotations_to_hash

  Arg 1      : hashref $vf_hash
  Arg 2      : string $custom_name
  Arg 3      : listref $annotations
  Example    : $hashref = $of->_add_custom_annotations_to_hash($vf, $vf_hash, $ex);
  Description: Adds custom annotation data to hash
  Returntype : hashref
  Exceptions : none
  Caller     : VariationFeatureOverlapAllele_to_output_hash(),
               VariationFeature_to_output_hash()
  Status     : Stable

=cut

sub _add_custom_annotations_to_hash {
  my ($self, $hash, $custom_name, $annots) = @_;

  foreach my $annot(@$annots) {
    push @{$hash->{$custom_name}}, $annot->{name};

    foreach my $field(keys %{$annot->{fields} || {}}) {
      push @{$hash->{$custom_name.'_'.$field}}, $annot->{fields}->{$field};
    }
  }

  return $hash;
}


=head2 VariationFeatureOverlapAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele $vfoa
  Arg 2      : hashref $vf_hash
  Arg 3      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : $hashref = $of->VariationFeatureOverlapAllele_to_output_hash($vfoa, $vf_hash, $vf);
  Description: Adds generic information to hash applicable to all
               VariationFeatureOverlapAllele sub-classes.
  Returntype : hashref
  Exceptions : none
  Caller     : TranscriptVariationAllele_to_output_hash(),
               RegulatoryFeatureVariationAllele_to_output_hash(),
               MotifFeatureVariationAllele_to_output_hash(),
               IntergenicVariationAllele_to_output_hash()
  Status     : Stable

=cut

sub VariationFeatureOverlapAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash, $vf) = @_;

  my @ocs = sort {$a->rank <=> $b->rank} @{$vfoa->get_all_OverlapConsequences};

  # consequence type(s)
  my $term_method = $self->{terms}.'_term';
  $hash->{Consequence} = [map {$_->$term_method || $_->SO_term} @ocs];

  # impact
  $hash->{IMPACT} = $ocs[0]->impact() if @ocs;

  # allele
  $hash->{Allele} = $vfoa->variation_feature_seq; 

  if(defined($vfoa->{shift_hash})&& defined($vfoa->{shift_hash}->{hgvs_allele_string}) && $self->param('shift_genomic'))
  {
    $hash->{Allele} = $vfoa->{shift_hash}->{hgvs_allele_string};
  }
  # allele number
  $hash->{ALLELE_NUM} = $vfoa->allele_number if $self->{allele_number};

  # reference allele
  $hash->{REF_ALLELE} = $vf->ref_allele_string if $self->{show_ref_allele};

  # picked?
  $hash->{PICK} = 1 if defined($vfoa->{PICK});

  # hgvs g.
  if($self->{hgvsg}) {
    $vf->{_hgvs_genomic} ||= $vf->hgvs_genomic($vf->slice, $self->{hgvsg_use_accession} ? undef : $vf->{chr});

    if(my $hgvsg = $vf->{_hgvs_genomic}->{$vfoa->variation_feature_seq}) {
      $hash->{HGVSg} = $hgvsg; 
    }
  }
  # spdi
  if($self->{spdi}) { 
    $vf->{_spdi_genomic} = $vf->spdi_genomic(); 
      
    if(my $spdi = $vf->{_spdi_genomic}->{$hash->{Allele}}){
      $hash->{SPDI} = $spdi;  
    }
  }

  # custom annotations
  foreach my $custom_name(keys %{$vf->{_custom_annotations} || {}}) {
    $self->_add_custom_annotations_to_hash(
      $hash,
      $custom_name,
      [
        grep {$_->{allele} && ($_->{allele} eq $hash->{Allele})}
        @{$vf->{_custom_annotations}->{$custom_name}}
      ]
    );
  }

  # colocated
  $self->add_colocated_variant_info($vf, $hash);

  # frequency data
  foreach my $ex(@{$vf->{existing} || []}) {
    $self->add_colocated_frequency_data($vf, $hash, $ex, $vfoa->{shift_hash});
  }

  return $hash;
}


=head2 BaseTranscriptVariationAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele $btva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->BaseTranscriptVariationAllele_to_output_hash($btva, $vf_hash);
  Description: Adds information to hash applicable to transcript.
  Returntype : hashref
  Exceptions : none
  Caller     : TranscriptVariationAllele_to_output_hash(),
               TranscriptStructuralVariationAllele_to_output_hash()
  Status     : Stable

=cut

sub BaseTranscriptVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  my $tv = $vfoa->base_variation_feature_overlap;
  my $tr = $tv->transcript;
  my $pre = $vfoa->_pre_consequence_predicates;

  # basics
  $hash->{Feature_type} = 'Transcript';
  $hash->{Feature}      = $tr->stable_id if $tr;
  $hash->{Feature}     .= '.'.$tr->version if $hash->{Feature} && $self->{transcript_version} && $hash->{Feature} !~ /\.\d+$/;

  # get gene
  $hash->{Gene} = $tr->{_gene_stable_id};

  # strand
  $hash->{STRAND} = $tr->strand + 0;

  my @attribs = @{$tr->get_all_Attributes()};

  # flags
  my @flags = grep {substr($_, 0, 4) eq 'cds_'} map {$_->{code}} @attribs;
  $hash->{FLAGS} = \@flags if scalar @flags;

  # exon/intron numbers
  if($self->{numbers}) {
    if($pre->{exon}) {
      if(my $num = $tv->exon_number) {
        $hash->{EXON} = $num;
      }
    }
    if($pre->{intron}) {
      if(my $num = $tv->intron_number) {
        $hash->{INTRON} = $num;
      }
    }
  }

  # protein domains
  if($self->{domains} && $pre->{coding}) {
    my $feats = $tv->get_overlapping_ProteinFeatures;

    my @strings;

    for my $feat (@$feats) {

      # do a join/grep in case self missing data
      my $label = join(':', grep {$_} ($feat->analysis->display_label, $feat->hseqname));

      # replace any special characters
      $label =~ s/[\s;=]/_/g;

      push @strings, $label;
    }

    $hash->{DOMAINS} = \@strings if @strings;
  }

  # distance to transcript
  if(grep {$DISTANCE_CONS{$_}} @{$hash->{Consequence} || []}) {
    $hash->{DISTANCE} = $tv->distance_to_transcript;
  }

  # gene symbol
  if($self->{symbol}) {
    my $symbol  = $tr->{_gene_symbol} || $tr->{_gene_hgnc};
    my $source  = $tr->{_gene_symbol_source};
    my $hgnc_id = $tr->{_gene_hgnc_id} if defined($tr->{_gene_hgnc_id});

    $hash->{SYMBOL} = $symbol if defined($symbol) && $symbol ne '-';
    ## encode spaces in symbol name for VCF
    $hash->{SYMBOL} =~ s/\s/\%20/g if $hash->{SYMBOL} && $self->{output_format} eq 'vcf' && !$self->{no_escape};

    $hash->{SYMBOL_SOURCE} = $source if defined($source) && $source ne '-';
    $hash->{HGNC_ID} = $hgnc_id if defined($hgnc_id) && $hgnc_id ne '-';
  }

  # CCDS
  $hash->{CCDS} = $tr->{_ccds} if
    $self->{ccds} &&
    defined($tr->{_ccds}) &&
    $tr->{_ccds} ne '-';

  # refseq xref
  $hash->{RefSeq} = [split(',', $tr->{_refseq})] if
    $self->{xref_refseq} &&
    defined($tr->{_refseq}) &&
    $tr->{_refseq} ne '-';

  # refseq match info
  if($self->{refseq} || $self->{merged}) {
    my @rseq_attrs = grep {$_->code =~ /^rseq/} @attribs;
    $hash->{REFSEQ_MATCH} = [map {$_->code} @rseq_attrs] if scalar @rseq_attrs;
  }

  if(my $status = $tr->{_bam_edit_status}) {
    $hash->{BAM_EDIT} = uc($status);
  }

  # protein ID
  $hash->{ENSP} = $tr->{_protein} if
    $self->{protein} &&
    defined($tr->{_protein}) &&
    $tr->{_protein} ne '-';

  # uniprot
  if($self->{uniprot}) {
    for my $db(qw(swissprot trembl uniparc uniprot_isoform)) {
      my $id = $tr->{'_'.$db};
      $id = undef if defined($id) && $id eq '-';
      $hash->{uc($db)} = [split(',', $id)] if defined($id);
    }
  }

  # canonical transcript
  $hash->{CANONICAL} = 'YES' if $self->{canonical} && $tr->is_canonical;

  # biotype
  $hash->{BIOTYPE} = $tr->biotype if $self->{biotype} && $tr->biotype;

  # source cache self transcript if using --merged
  $hash->{SOURCE} = $tr->{_source_cache} if defined $tr->{_source_cache};

  # gene phenotype
  $hash->{GENE_PHENO} = 1 if $self->{gene_phenotype} && $tr->{_gene_phenotype};
  if($self->{mane_select} && (my ($mane) = grep {$_->code eq 'MANE_Select'} @attribs)) {
    if(my $mane_value = $mane->value) {
      $hash->{MANE_SELECT} = $mane_value;
    }
  }

  if($self->{mane_plus_clinical} && (my ($mane) = grep {$_->code eq 'MANE_Plus_Clinical'} @attribs)) {
    if(my $mane_value = $mane->value) {
      $hash->{MANE_PLUS_CLINICAL} = $mane_value;
    }
  }
  
  # transcript support level
  if($self->{tsl} && (my ($tsl) = grep {$_->code eq 'TSL'} @attribs)) {
    if($tsl->value =~ m/tsl(\d+)/) {
      $hash->{TSL} = $1 if $1;
    }
  }

  # APPRIS
  if($self->{appris} && (my ($appris) = grep {$_->code eq 'appris'} @attribs)) {
    if(my $value = $appris->value) {
      $value =~ s/principal/P/;
      $value =~ s/alternative/A/;
      $hash->{APPRIS} = $value;
    }
  }

  # miRNA structure
  if($self->{mirna} && (my ($mirna_attrib) = grep {$_->code eq 'ncRNA'} @attribs)) {
    my ($start, $end, $struct) = split /\s+|\:/, $mirna_attrib->value;

    my ($cdna_start, $cdna_end) = ($tv->cdna_start, $tv->cdna_end);

    if(
      defined($struct) && $struct =~ /[\(\.\)]+/ &&
      $start && $end && $cdna_start && $cdna_end &&
      overlap($start, $end, $cdna_start, $cdna_end)
    ) {

      # account for insertions
      ($cdna_start, $cdna_end) = ($cdna_end, $cdna_start) if $cdna_start > $cdna_end;
    
      # parse out structure
      my @struct;
      while($struct =~ m/([\.\(\)])([0-9]+)?/g) {
        my $num = $2 || 1;
        push @struct, $1 for(1..$num);
      }
      
      # get struct element types overlapped by variant
      my %chars;
      for my $pos($cdna_start..$cdna_end) {
        $pos -= $start;
        next if $pos < 0 or $pos > scalar @struct;
        $chars{$struct[$pos]} = 1;
      }
      
      # map element types to SO terms
      my %map = (
        '(' => 'miRNA_stem',
        ')' => 'miRNA_stem',
        '.' => 'miRNA_loop'
      );
      
      $hash->{miRNA} = [sort map {$map{$_}} keys %chars];
    }
  }

  return $hash;
}


=head2 TranscriptVariationAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::TranscriptVariationAllele $tva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->TranscriptVariationAllele_to_output_hash($tva, $vf_hash);
  Description: Adds information to hash applicable to transcript and allele combination.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub TranscriptVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  # run "super" methods
  $hash = $self->VariationFeatureOverlapAllele_to_output_hash(@_);
  $hash = $self->BaseTranscriptVariationAllele_to_output_hash(@_);
  
  my $shift_length = (defined($vfoa->{shift_hash}) ? $vfoa->{shift_hash}->{shift_length} : 0);
  $shift_length ||= 0;
  
  return undef unless $hash;

  $hash->{SHIFT_LENGTH} = $shift_length if ($self->param('shift_3prime') && $self->param('shift_length')); 

  my $tv = $vfoa->base_variation_feature_overlap;
  my $tr = $tv->transcript;
  
  my $strand = defined($tr->strand) ? $tr->strand : 1;

  if($self->{shift_genomic})
  {
    my $vf = $vfoa->variation_feature;
    $hash->{Location} = ($vf->{chr} || $vf->seq_region_name).':'.format_coords($vf->{start} + ($shift_length * $strand), $vf->{end} + ($shift_length * $strand));
  }

  my $vep_cache = $tr->{_variation_effect_feature_cache};

  my $pre = $vfoa->_pre_consequence_predicates();

  if($pre->{within_feature}) {

    # exonic only
    if($pre->{exon}) {
      my $shifting_offset = $shift_length * $strand;
      $hash->{cDNA_position}  = format_coords($tv->cdna_start(undef,$shifting_offset), $tv->cdna_end(undef,$shifting_offset));  

      $hash->{cDNA_position} .= '/'.$tr->length if $self->{total_length};

      # coding only
      if($pre->{coding}) {

        $hash->{Amino_acids} = $vfoa->pep_allele_string;
        $hash->{Codons}      = $vfoa->display_codon_allele_string;
        $shifting_offset = 0 if defined($tv->{_boundary_shift}) && $tv->{_boundary_shift} == 1;

        $hash->{CDS_position}  = format_coords($tv->cds_start, $tv->cds_end);
        $hash->{CDS_position} .= '/'.length($vep_cache->{translateable_seq})
          if $self->{total_length} && $vep_cache->{translateable_seq};

        $hash->{Protein_position}  = format_coords($tv->translation_start(undef, $shifting_offset), $tv->translation_end(undef, $shifting_offset));
        $hash->{Protein_position} .= '/'.length($vep_cache->{peptide})
          if $self->{total_length} && $vep_cache->{peptide};

        $self->add_sift_polyphen($vfoa, $hash);
      }
    }
    my $strand = $tr->strand() > 0 ? 1 : -1;
    # HGVS
    if($self->{hgvsc}) {
      my $hgvs_t = $vfoa->hgvs_transcript(undef, !$self->param('shift_3prime'));
      my $offset = defined($vfoa->{shift_hash}) ? $vfoa->{shift_hash}->{_hgvs_offset} : 0;

      $hash->{HGVSc} = $hgvs_t if $hgvs_t;
      $hash->{HGVS_OFFSET} = $offset * $strand if $offset && $hgvs_t;
    }

    if($self->{hgvsp}) {
      $vfoa->{remove_hgvsp_version} = 1 if $self->{remove_hgvsp_version};
      my $hgvs_p = $vfoa->hgvs_protein;
      my $offset = $vfoa->hgvs_offset;

      # URI encode "="
      $hgvs_p =~ s/\=/\%3D/g if $hgvs_p && !$self->{no_escape};

      $hash->{HGVSp} = $hgvs_p if $hgvs_p;
      $hash->{HGVS_OFFSET} = $offset * $strand if $offset && $hgvs_p;
    }
    $hash->{REFSEQ_OFFSET} = $vfoa->{refseq_misalignment_offset} if defined($vfoa->{refseq_misalignment_offset}) && $vfoa->{refseq_misalignment_offset} != 0;
  }



  if($self->{use_transcript_ref}) {
    my $ref_tva = $tv->get_reference_TranscriptVariationAllele;
    $hash->{USED_REF} = $ref_tva->variation_feature_seq;
    $hash->{USED_REF} = $ref_tva->{shift_hash}->{ref_orig_allele_string} if !$self->{shift_3prime} && defined($ref_tva->{shift_hash});
    $hash->{GIVEN_REF} = $ref_tva->{given_ref};
  }

  return $hash;
}


=head2 add_sift_polyphen

  Arg 1      : Bio::EnsEMBL::Variation::TranscriptVariationAllele $tva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->add_sift_polyphen($tva, $vf_hash);
  Description: Adds SIFT and PolyPhen data to hash.
  Returntype : hashref
  Exceptions : none
  Caller     : TranscriptVariationAllele_to_output_hash()
  Status     : Stable

=cut

sub add_sift_polyphen {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  foreach my $tool (qw(SIFT PolyPhen)) {
    my $lc_tool = lc($tool);

    if (my $opt = $self->{$lc_tool}) {
      my $want_pred  = $opt =~ /^p/i;
      my $want_score = $opt =~ /^s/i;
      my $want_both  = $opt =~ /^b/i;

      if ($want_both) {
        $want_pred  = 1;
        $want_score = 1;
      }

      next unless $want_pred || $want_score;

      my $pred_meth  = $lc_tool.'_prediction';
      my $score_meth = $lc_tool.'_score';
      my $analysis   = $self->{polyphen_analysis} if $lc_tool eq 'polyphen';

      my $pred = $vfoa->$pred_meth($analysis);

      if($pred) {

        if ($want_pred) {
          $pred =~ s/\s+/\_/g;
          $pred =~ s/\_\-\_/\_/g;
          $hash->{$tool} = $pred;
        }

        if ($want_score) {
          my $score = $vfoa->$score_meth($analysis);

          if(defined $score) {
            if($want_pred) {
              $hash->{$tool} .= "($score)";
            }
            else {
              $hash->{$tool} = $score;
            }
          }
        }
      }

      # update stats
      $self->stats->log_sift_polyphen($tool, $pred) if $pred && !$self->{no_stats};
    }
  }

  return $hash;
}


=head2 RegulatoryFeatureVariationAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele $rfva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->RegulatoryFeatureVariationAllele_to_output_hash($rfva, $vf_hash);
  Description: Adds information to hash applicable to regulatory feature and allele combination.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub RegulatoryFeatureVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  # run "super" method
  $hash = $self->VariationFeatureOverlapAllele_to_output_hash(@_);
  return undef unless $hash;

  my $rf = $vfoa->regulatory_feature;

  $hash->{Feature_type} = 'RegulatoryFeature';
  $hash->{Feature}      = $rf->stable_id;
  $hash->{CELL_TYPE}    = $self->get_cell_types($rf) if $self->{cell_type};

  $hash->{BIOTYPE} = ref($rf->{feature_type}) ? $rf->{feature_type}->{so_name} : $rf->{feature_type} if defined($rf->{feature_type});

  return $hash;
}


=head2 MotifFeatureVariationAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::MotifFeatureVariationAllele $mfva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->MotifFeatureVariationAllele_to_output_hash($mfva, $vf_hash);
  Description: Adds information to hash applicable to motif feature and allele combination.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub MotifFeatureVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  # run "super" method
  $hash = $self->VariationFeatureOverlapAllele_to_output_hash(@_);
  return undef unless $hash;

  my $mf = $vfoa->motif_feature;

  # check that the motif has a binding matrix, if not there's not
  # much we can do so don't return anything
  return undef unless defined $mf->get_BindingMatrix;
  my $matrix = $mf->get_BindingMatrix;
  my $matrix_id = $matrix->stable_id;

  my $mf_stable_id = $mf->stable_id;
  $hash->{Feature_type} = 'MotifFeature';
  $hash->{Feature}      = $mf_stable_id;
  $hash->{MOTIF_NAME}   = $matrix_id;
  my @transcription_factors = ();
  my $associated_transcription_factor_complexes = $matrix->{associated_transcription_factor_complexes};
  foreach my $tfc (@{$associated_transcription_factor_complexes}) {
    push @transcription_factors, $tfc->{display_name};
  }
  $hash->{TRANSCRIPTION_FACTORS} = \@transcription_factors;
  $hash->{STRAND}       = $mf->strand + 0;
  $hash->{CELL_TYPE}    = $self->get_cell_types($mf) if $self->{cell_type};
  $hash->{MOTIF_POS}    = $vfoa->motif_start if defined $vfoa->motif_start;
  $hash->{HIGH_INF_POS} = ($vfoa->in_informative_position ? 'Y' : 'N');

  my $delta = $vfoa->motif_score_delta if $vfoa->variation_feature_seq =~ /^[ACGT]+$/;
  $hash->{MOTIF_SCORE_CHANGE} = sprintf("%.3f", $delta) if defined $delta;

  return $hash;
}


=head2 get_cell_types

  Arg 1      : Bio::EnsEMBL::Funcgen::MotifFeature or Bio::EnsEMBL::Funcgen::RegulatoryFeature $ft
  Example    : $cell_types = $of->get_cell_types($ft);
  Description: Gets cell types for this regulatory or motif feature
  Returntype : listref of strings
  Exceptions : none
  Caller     : MotifFeatureVariationAllele_to_output_hash(),
               RegulatoryFeatureVariationAllele_to_output_hash()
  Status     : Stable

=cut

sub get_cell_types {
  my $self = shift;
  my $ft = shift;

  return [
    map {s/\s+/\_/g; $_}
    map {$_.':'.$ft->{cell_types}->{$_}}
    grep {$ft->{cell_types}->{$_}}
    @{$self->{cell_type}}
  ];
}


=head2 IntergenicVariationAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::IntergenicVariationAllele $iva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->IntergenicVariationAllele_to_output_hash($iva, $vf_hash);
  Description: Just a placeholder really as no extra information is added for
               intergenic variants.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub IntergenicVariationAllele_to_output_hash {
  my $self = shift;
  my $iva = shift;
  my $hash = shift;
    
  ## By default, shifting in the 3' direction will not update the 'location' field unless '--shift_genomic' is supplied
  unless(!$self->{shift_3prime} || !$self->param('shift_genomic')) { #if shifting without shift genomic, we still want to run $iva->genomic_shift
    $iva->genomic_shift;
    if (defined($iva->{shift_hash}->{shift_length}) && $self->param('shift_genomic')) {
      my $vf = $iva->variation_feature;
      $hash->{Location} = ($vf->{chr} || $vf->seq_region_name).':'.
        format_coords($vf->{start} + $iva->{shift_hash}->{shift_length}, $vf->{end} + $iva->{shift_hash}->{shift_length});
    }
  }
  
  return $self->VariationFeatureOverlapAllele_to_output_hash($iva, $hash, @_);
}


### SV-type methods
###################


=head2 BaseStructuralVariationOverlapAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::BaseStructuralVariationOverlapAllele $bsvoa
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->BaseStructuralVariationOverlapAllele_to_output_hash($bsvoa, $vf_hash);
  Description: Adds basic information to hash applicable to all
               StructuralVariationOverlapAllele classes.
  Returntype : hashref
  Exceptions : none
  Caller     : StructuralVariationOverlapAllele_to_output_hash()
  Status     : Stable

=cut

sub BaseStructuralVariationOverlapAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  my $svf = $vfoa->base_variation_feature;

  $hash->{Allele} = $svf->class_SO_term;

  # allele number
  $hash->{ALLELE_NUM} = $vfoa->allele_number if $self->{allele_number};

  my @ocs = sort {$a->rank <=> $b->rank} @{$vfoa->get_all_OverlapConsequences};

  # consequence type(s)
  my $term_method = $self->{terms}.'_term';
  $hash->{Consequence} = [map {$_->$term_method || $_->SO_term} @ocs];

  # impact
  $hash->{IMPACT} = $ocs[0]->impact() if @ocs;

  # allele number
  # if($self->{allele_number}) {
  #   $hash->{ALLELE_NUM} = $vfoa->allele_number if $vfoa->can('allele_number');
  # }

  # picked?
  $hash->{PICK} = 1 if defined($vfoa->{PICK});

  return $hash;
}


=head2 StructuralVariationOverlapAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele $svoa
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->StructuralVariationOverlapAllele_to_output_hash($svoa, $vf_hash);
  Description: Adds basic information to hash applicable to all
               StructuralVariationOverlapAllele classes.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub StructuralVariationOverlapAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  # run "super" method
  $hash = $self->BaseStructuralVariationOverlapAllele_to_output_hash(@_);
  return undef unless $hash;

  my $feature = $vfoa->feature;
  my $svf = $vfoa->base_variation_feature;

  # get feature type
  my $feature_type = (split '::', ref($feature))[-1];

  $hash->{Feature_type} = $feature_type;
  $hash->{Feature}      = $feature->stable_id;

  # work out overlap amounts
  my $overlap_start  = (sort {$a <=> $b} ($svf->start, $feature->start))[-1];
  my $overlap_end    = (sort {$a <=> $b} ($svf->end, $feature->end))[0];
  my $overlap_length = ($overlap_end - $overlap_start) + 1;
  my $overlap_pc     = 100 * ($overlap_length / (($feature->end - $feature->start) + 1));

  $hash->{OverlapBP} = $overlap_length if $overlap_length > 0;
  $hash->{OverlapPC} = sprintf("%.2f", $overlap_pc) if $overlap_pc > 0;

  # cell types
  $hash->{CELL_TYPE} = $self->get_cell_types($feature)
    if $self->{cell_type} && $feature_type =~ /(Motif|Regulatory)Feature/;

  return $hash;
}


=head2 TranscriptStructuralVariationAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele $tsva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->TranscriptStructuralVariationAllele_to_output_hash($tsva, $vf_hash);
  Description: Adds transcript information to hash for SV overlaps.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub TranscriptStructuralVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  # run "super" methods
  $hash = $self->StructuralVariationOverlapAllele_to_output_hash(@_);
  $hash = $self->BaseTranscriptVariationAllele_to_output_hash(@_);
  return undef unless $hash;

  my $svo = $vfoa->base_variation_feature_overlap;
  my $pre = $vfoa->_pre_consequence_predicates;
  my $tr  = $svo->transcript;
  my $vep_cache = $tr->{_variation_effect_feature_cache};

  if($pre->{within_feature}) {

    # exonic only
    if($pre->{exon}) {

      $hash->{cDNA_position}  = format_coords($svo->cdna_start, $svo->cdna_end);
      $hash->{cDNA_position} .= '/'.$tr->length if $self->{total_length};

      # coding only
      if($pre->{coding}) {

        $hash->{CDS_position}  = format_coords($svo->cds_start, $svo->cds_end);
        $hash->{CDS_position} .= '/'.length($vep_cache->{translateable_seq})
          if $self->{total_length} && $vep_cache->{translateable_seq};

        $hash->{Protein_position}  = format_coords($svo->translation_start, $svo->translation_end);
        $hash->{Protein_position} .= '/'.length($vep_cache->{peptide})
          if $self->{total_length} && $vep_cache->{peptide};
      }
    }
  }

  return $hash;
}


=head2 IntergenicStructuralVariationAllele_to_output_hash

  Arg 1      : Bio::EnsEMBL::Variation::IntergenicStructuralVariationAllele $isva
  Arg 2      : hashref $vf_hash
  Example    : $hashref = $of->IntergenicStructuralVariationAllele_to_output_hash($isva, $vf_hash);
  Description: Just a placeholder really as no extra information is added for
               intergenic variants.
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub IntergenicStructuralVariationAllele_to_output_hash {
  my $self = shift;
  return $self->BaseStructuralVariationOverlapAllele_to_output_hash(@_);
}


### Other methods
#################


=head2 plugins

  Example    : $plugins = $of->plugins();
  Description: Gets all plugins to be run.
  Returntype : listref of Bio::EnsEMBL::Variation::Utils::BaseVepPlugin
  Exceptions : none
  Caller     : run_plugins(), get_plugin_headers()
  Status     : Stable

=cut

sub plugins {
  return $_[0]->{plugins} || [];
}


=head2 run_plugins

  Arg 1      : Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele $bvfoa
  Arg 2      : hashref $vf_hash
  Arg 3      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : $vf_hash = $of->run_plugins($bvfoa, $vf_hash, $vf);
  Description: Runs all plugins for given BaseVariationFeatureOverlapAllele
  Returntype : hashref
  Exceptions : none
  Caller     : get_all_VariationFeatureOverlapAllele_output_hashes()
  Status     : Stable

=cut

sub run_plugins {
  my ($self, $bvfoa, $line_hash, $vf) = @_;

  my $skip_line = 0;

  for my $plugin (@{ $self->plugins || [] }) {

    # check that this plugin is interested in this type of variation feature
    if($plugin->check_variant_feature_type(ref($vf || $bvfoa->base_variation_feature))) {

      # check that this plugin is interested in this type of feature
      if($plugin->check_feature_type(ref($bvfoa->feature) || 'Intergenic')) {

        eval {
          my $plugin_results = $plugin->run($bvfoa, $line_hash);

          if (defined $plugin_results) {
            if (ref $plugin_results eq 'HASH') {
              for my $key (keys %$plugin_results) {
                $line_hash->{$key} = $plugin_results->{$key};
              }
            }
            else {
              $self->warning_msg("Plugin '".(ref $plugin)."' did not return a hashref, output ignored!");
            }
          }
          else {
            # if a plugin returns undef, that means it want to filter out this line
            $skip_line = 1;
          }
        };
        if($@) {
          $self->warning_msg("Plugin '".(ref $plugin)."' went wrong: $@");
        }

        # there's no point running any other plugins if we're filtering this line,
        # because the first plugin to skip the line wins, so we might as well last
        # out of the loop now and avoid any unnecessary computation
        last if $skip_line;
      }
    }
  }

  return $skip_line ? undef : $line_hash;
}


=head2 rejoin_variants_in_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $of->rejoin_variants_in_InputBuffer($ib);
  Description: Rejoins linked variants that were split in InputBuffer, typically
               by user providing --minimal flag.
  Returntype : none
  Exceptions : none
  Caller     : get_all_output_hashes_by_InputBuffer()
  Status     : Stable

=cut

sub rejoin_variants_in_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my @joined_list = ();

  # backup stat logging status
  my $no_stats = $self->{no_stats};
  $self->{no_stats} = 1;

  foreach my $vf(@{$buffer->buffer}) {

    # reset original one
    if(defined($vf->{original_allele_string})) {

      # do consequence stuff
      $self->get_all_output_hashes_by_VariationFeature($vf);

      $vf->{allele_string}    = $vf->{original_allele_string};
      $vf->{seq_region_start} = $vf->{start} = $vf->{original_start};
      $vf->{seq_region_end}   = $vf->{end}   = $vf->{original_end};

      push @joined_list, $vf;
    }

    # this one needs to be merged in
    elsif(defined($vf->{merge_with})) {
      my $original = $vf->{merge_with};

      # do consequence stuff
      $self->get_all_output_hashes_by_VariationFeature($vf);

      # now we have to copy the [Feature]Variation objects
      # we can't simply copy the alleles as the coords will be different
      # better to make new keys
      # we also have to set the VF pointer to the original

      # copy transcript variations etc
      foreach my $type(map {$_.'_variations'} qw(transcript motif_feature regulatory_feature)) {
        foreach my $key(keys %{$vf->{$type} || {}}) {
          my $val = $vf->{$type}->{$key};
          $val->base_variation_feature($original);
          # rename the key they're stored under
          $original->{$type}->{$vf->{allele_string}.'_'.$key} = $val;
        }
      }

      # intergenic variation is a bit different
      # there is only one, and no reference feature to key on
      # means we have to copy over alleles manually
      if(my $iv = $vf->{intergenic_variation}) {
        my $bfvo = $original->{intergenic_variation};
        $iv->base_variation_feature($original);

        if(my $oiv = $original->{intergenic_variation}) {
          foreach my $alt (@{$iv->{alt_alleles}}) {
            $alt->base_variation_feature_overlap($bfvo);
            push @{$oiv->{alt_alleles}}, $alt;
          }
          $oiv->{_alleles_by_seq}->{$_->variation_feature_seq} = $_ for @{$oiv->{alt_alleles}};
        }

        # this probably won't happen, but can't hurt to cover all bases
        else {
            $original->{intergenic_variation} = $iv;
        }
      }

      # reset these keys, they can be recalculated
      delete $original->{$_} for qw(overlap_consequences _most_severe_consequence);
    }

    # normal
    else {
      push @joined_list, $vf;
    }
  }

  $self->{no_stats} = $no_stats;

  $buffer->buffer(\@joined_list);
}


### Header etc methods
### These will be called when generating the actual output
##########################################################


=head2 header_info

  Example    : $info = $of->header_info();
  Description: Retrieve header info passed in on object creation
  Returntype : hashref
  Exceptions : none
  Caller     : get_custom_headers(), sub-classes
  Status     : Stable

=cut

sub header_info {
  my $self = shift;
  return $self->{header_info} || {};
}


=head2 headers

  Example    : $headers = $of->headers();
  Description: Get list of headers to print out for this output format.
               This is a stub method, overridden in child classes.
  Returntype : listref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub headers {
  return [];
}


=head2 get_plugin_headers

  Example    : $headers = $of->get_plugin_headers();
  Description: Get headers from plugins
  Returntype : arrayref of arrayrefs [$key, $header]
  Exceptions : none
  Caller     : description_headers() in child classes
  Status     : Stable

=cut

sub get_plugin_headers {
  my $self = shift;

  my @headers = ();

  for my $plugin (@{$self->plugins}) {
    if (my $hdr = $plugin->get_header_info) {
      for my $key (sort keys %$hdr) {
        push @headers, [$key, $hdr->{$key}];
      }
    }
  }

  return \@headers;
}


=head2 get_custom_headers

  Example    : $headers = $of->get_custom_headers();
  Description: Get headers from custom data files
  Returntype : arrayref of arrayrefs [$key, $header]
  Exceptions : none
  Caller     : description_headers() in child classes
  Status     : Stable

=cut

sub get_custom_headers {
  my $self = shift;

  my @headers;

  foreach my $custom(@{$self->header_info->{custom_info} || []}) {
    push @headers, [$custom->{short_name}, sprintf("%s (%s)", $custom->{file}, $custom->{type})];
    
    foreach my $field(@{$custom->{fields} || []}) {
      push @headers, [
        sprintf("%s_%s", $custom->{short_name}, $field),
        sprintf("%s field from %s", $field, $custom->{file})
      ];
    }
  }

  return \@headers;
}


=head2 flag_fields

  Example    : $fields = $of->flag_fields();
  Description: Get list of fields that will appear in output given user config
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : fields() in child classes
  Status     : Stable

=cut

sub flag_fields {
  my $self = shift;
  
  # get all fields
  my @tmp = (
    map {@{$_->{fields}}}
    map {$_->[0]}
    grep {
      ref($_->[1]) eq 'ARRAY' ? scalar @{$_->[1]} : $_->[1]
    }
    map {[$_, $self->param($_->{flag})]}
    @Bio::EnsEMBL::VEP::Constants::FLAG_FIELDS
  );

  # uniquify, retaining order
  my (%seen, @return);
  for(@tmp) {
    push @return, $_ unless $seen{$_};
    $seen{$_} = 1;
  }

  return \@return;
}

1;
