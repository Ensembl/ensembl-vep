=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::OutputFactory;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Scalar::Util qw(looks_like_number);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::Constants;
use Bio::EnsEMBL::VEP::Utils qw(format_coords);
use Bio::EnsEMBL::VEP::Constants;

use Bio::EnsEMBL::VEP::OutputFactory::VEP_output;
use Bio::EnsEMBL::VEP::OutputFactory::VCF;
use Bio::EnsEMBL::VEP::OutputFactory::JSON;
use Bio::EnsEMBL::VEP::OutputFactory::Tab;

my %SO_RANKS = map {$_->SO_term => $_->rank} values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;

my %FORMAT_MAP = (
  'vcf'     => 'VCF',
  'ensembl' => 'VEP_output',
  'vep'     => 'VEP_output',
  'tab'     => 'Tab',
  'json'    => 'JSON',
);

my %DISTANCE_CONS = (upstream_gene_variant => 1, downstream_gene_variant => 1);

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
    no_escape
    pick_order
    allele_number

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
    gmaf
    maf_1kg
    maf_esp
    maf_exac
    pubmed

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
    tsl
    appris
    gene_phenotype

    total_length
    hgvs
    sift
    polyphen
    polyphen_analysis
  )]);

  my $hashref = $_[0];

  if(my $format = $hashref->{format}) {

    delete $hashref->{format};

    $format = lc($format);
    throw("ERROR: Unknown or unsupported output format $format\n") unless $FORMAT_MAP{$format};

    my $class = 'Bio::EnsEMBL::VEP::OutputFactory::'.$FORMAT_MAP{$format};
    return $class->new({%$hashref, config => $self->config});
  }

  $self->{header_info} = $hashref->{header_info} if $hashref->{header_info};

  $self->{plugins} = $hashref->{plugins} if $hashref->{plugins};

  return $self;
}

# sub print_InputBuffer {
#   my $self = shift;
#   my $buffer = shift;

#   $self->print_line($_) for @{$self->get_all_lines_by_InputBuffer($buffer)};
# }

# this is the method that will be called by Runner
# may be re-implemented by child classes to account for differing structure
sub get_all_lines_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;
  return [map {$self->output_hash_to_line($_)} @{$self->get_all_output_hashes_by_InputBuffer($buffer)}];
}

# this method can also be called by Runner if unprocessed output is required
sub get_all_output_hashes_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  $self->rejoin_variants_in_InputBuffer($buffer) if $buffer->rejoin_required;

  return [
    map {@{$self->get_all_output_hashes_by_VariationFeature($_)}}
    @{$buffer->buffer}
  ];
}

# method to retrieve header info passed in on object creation
sub header_info {
  my $self = shift;
  return $self->{header_info} || {};
}

# sub-classes will typically re-implement this
sub headers {
  return [];
}

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

# this method returns all the fields that will be populated given user parameters
sub flag_fields {
  my $self = shift;
  return [
    map {@{$_->{fields}}}
    map {$_->[0]}
    grep {
      ref($_->[1]) eq 'ARRAY' ? scalar @{$_->[1]} : $_->[1]
    }
    map {[$_, $self->param($_->{flag})]}
    @Bio::EnsEMBL::VEP::Constants::FLAG_FIELDS
  ];
}

# this wrapper method does the following:
# 1) fetches all the VariationFeatureOverlapAllele objects
# 2) filters them if any filtering options are applied
# 3) converts them into a flat hashref with the relevant output data
sub get_all_output_hashes_by_VariationFeature {
  my $self = shift;
  my $vf = shift;

  # get the basic initial hash; this basically contains the location and ID
  my $hash = $self->VariationFeature_to_output_hash($vf);

  return $self->get_all_VariationFeatureOverlapAllele_output_hashes($vf, $hash);
}

# this method is separated out from get_all_output_hashes_by_VariationFeature
# as the JSON sub-class needs this data _without_ that from VariationFeature_to_output_hash merged in
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

sub get_all_VariationFeatureOverlapAlleles {
  my $self = shift;
  my $vf = shift;

  # no intergenic?
  return [] if $self->{no_intergenic} && defined($vf->{intergenic_variation});

  # get all VFOAs
  # need to be sensitive to whether --coding_only is switched on
  my $vfos;

  # if coding only just get transcript & intergenic ones
  if($self->{coding_only}) {
    @$vfos = grep {defined($_)} (
      @{$vf->get_all_TranscriptVariations},
      $vf->get_IntergenicVariation
    );
  }
  else {
    $vfos = $vf->get_all_VariationFeatureOverlaps;
  }

  # grep out non-coding?
  @$vfos = grep {$_->can('affects_cds') && $_->affects_cds} @$vfos if $self->{coding_only};

  # method name stub for getting *VariationAlleles
  my $allele_method = $self->{process_ref_homs} ? 'get_all_' : 'get_all_alternate_';  
  my $method = $allele_method.'VariationFeatureOverlapAlleles';

  return $self->filter_VariationFeatureOverlapAlleles([map {@{$_->$method}} @{$vfos}]);
}

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

  return $self->filter_VariationFeatureOverlapAlleles([map {@{$_->get_all_StructuralVariationOverlapAlleles}} @{$vfos}]);
}

# this method chooses the appropriate filtering method
# the "flag" options are not strictly filters, but they operate on the same logic
sub filter_VariationFeatureOverlapAlleles {
  my $self = shift;
  my $vfoas = shift;

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

# picks the worst self a list self VariationFeatureOverlapAlleles
# VFOAs are ordered by a heirarchy:
# 1: canonical
# 2: transcript support level
# 3: biotype (protein coding favoured)
# 4: consequence rank
# 5: transcript length
# 6: transcript from Ensembl?
# 7: transcript from RefSeq?
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

# pick one VariationFeatureOverlapAllele per gene
# allow non-transcript types to pass through
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

# initialize an output hash for a VariationFeature
sub VariationFeature_to_output_hash {
  my $self = shift;
  my $vf = shift;

  my $hash = {
    Uploaded_variation  => $vf->variation_name,
    Location            => ($vf->{chr} || $vf->seq_region_name).':'.format_coords($vf->{start}, $vf->{end}),
  };

  # overlapping variants
  $self->add_colocated_variant_info($vf, $hash);

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

  $self->stats->log_VariationFeature($vf, $hash) unless $self->{no_stats};

  return $hash;
}

sub add_colocated_variant_info {
  my $self = shift;
  my $vf = shift;
  my $hash = shift;

  return unless $vf->{existing} && scalar @{$vf->{existing}};

  my @existing = @{$vf->{existing}};

  my $tmp = {};

  foreach my $ex(@{$vf->{existing}}) {

    # ID
    push @{$tmp->{Existing_variation}}, $ex->{variation_name} if $ex->{variation_name};

    # GMAF?
    if($self->{gmaf}) {
      push @{$tmp->{GMAF}}, $ex->{minor_allele}.':'.$ex->{minor_allele_freq}
        if $ex->{minor_allele} && looks_like_number($ex->{minor_allele_freq});
    }

    # other freqs we can treat all the same
    my @pops = ();
    push @pops, qw(AFR AMR ASN EAS EUR SAS)                                    if $self->{maf_1kg};
    push @pops, qw(AA EA)                                                      if $self->{maf_esp};
    push @pops, ('ExAC', map {'ExAC_'.$_} qw(Adj AFR AMR EAS FIN NFE OTH SAS)) if $self->{maf_exac};

    push @{$tmp->{$_.'_MAF'}}, $ex->{$_} for grep {defined($ex->{$_})} @pops;

    # clin sig and pubmed?
    push @{$tmp->{CLIN_SIG}}, $ex->{clin_sig} if $ex->{clin_sig};
    push @{$tmp->{PUBMED}}, $ex->{pubmed} if $self->{pubmed} && $ex->{pubmed};

    # somatic?
    push @{$tmp->{SOMATIC}}, $ex->{somatic} ? 1 : 0;

    # phenotype or disease
    push @{$tmp->{PHENO}}, $ex->{phenotype_or_disease} ? 1 : 0;
  }

  # post-process to remove all-0, e.g. SOMATIC
  foreach my $key(keys %$tmp) {
    delete $tmp->{$key} unless grep {$_} @{$tmp->{$key}};
  }

  # copy to hash
  $hash->{$_} = $tmp->{$_} for keys %$tmp;

  # frequencies used to filter will appear here
  if($vf->{_freq_check_freqs}) {
    my @freqs;

    foreach my $p(keys %{$vf->{_freq_check_freqs}}) {
      foreach my $a(keys %{$vf->{_freq_check_freqs}->{$p}}) {
        push @freqs, sprintf(
          '%s:%s:%g',
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

sub _add_custom_annotations_to_hash {
  my ($self, $hash, $custom_name, $annots) = @_;

  foreach my $annot(@$annots) {
    $hash->{$custom_name} = $annot->{name};

    foreach my $field(keys %{$annot->{fields} || {}}) {
      $hash->{$custom_name.'_'.$field} = $annot->{fields}->{$field};
    }
  }
}

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

  # allele number
  $hash->{ALLELE_NUM} = $vfoa->allele_number if $self->{allele_number};

  # picked?
  $hash->{PICK} = 1 if defined($vfoa->{PICK});

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

  return $hash;
}

sub BaseTranscriptVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  my $tv = $vfoa->base_variation_feature_overlap;
  my $tr = $tv->transcript;
  my $pre = $vfoa->_pre_consequence_predicates;

  # basics
  $hash->{Feature_type} = 'Transcript';
  $hash->{Feature}      = $tr->stable_id if $tr;

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
    $hash->{SYMBOL_SOURCE} = $source if defined($source) && $source ne '-';
    $hash->{HGNC_ID} = $hgnc_id if defined($hgnc_id) && $hgnc_id ne '-';
  }

  # CCDS
  $hash->{CCDS} = $tr->{_ccds} if
    $self->{ccds} &&
    defined($tr->{_ccds}) &&
    $tr->{_ccds} ne '-';

  # refseq xref
  $hash->{RefSeq} = $tr->{_refseq} if
    $self->{xref_refseq} &&
    defined($tr->{_refseq}) &&
    $tr->{_refseq} ne '-';

  # refseq match info
  if($self->{refseq} || $self->{merged}) {
    my @rseq_attrs = grep {$_->code =~ /^rseq/} @attribs;
    $hash->{REFSEQ_MATCH} = [map {$_->code} @rseq_attrs] if scalar @rseq_attrs;
  }

  # protein ID
  $hash->{ENSP} = $tr->{_protein} if
    $self->{protein} &&
    defined($tr->{_protein}) &&
    $tr->{_protein} ne '-';

  # uniprot
  if($self->{uniprot}) {
    for my $db(qw(swissprot trembl uniparc)) {
      my $id = $tr->{'_'.$db};
      $id = undef if defined($id) && $id eq '-';
      $hash->{uc($db)} = $id if defined($id);
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

  return $hash;
}

sub TranscriptVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  # run "super" methods
  $hash = $self->VariationFeatureOverlapAllele_to_output_hash(@_);
  $hash = $self->BaseTranscriptVariationAllele_to_output_hash(@_);
  return undef unless $hash;

  my $tv = $vfoa->base_variation_feature_overlap;
  my $tr = $tv->transcript;
  my $vep_cache = $tr->{_variation_effect_feature_cache};

  my $pre = $vfoa->_pre_consequence_predicates();

  if($pre->{within_feature}) {

    # exonic only
    if($pre->{exon}) {

      $hash->{cDNA_position}  = format_coords($tv->cdna_start, $tv->cdna_end);
      $hash->{cDNA_position} .= '/'.$tr->length if $self->{total_length};

      # coding only
      if($pre->{coding}) {

        $hash->{Amino_acids} = $vfoa->pep_allele_string;
        $hash->{Codons}      = $vfoa->display_codon_allele_string;

        $hash->{CDS_position}  = format_coords($tv->cds_start, $tv->cds_end);
        $hash->{CDS_position} .= '/'.length($vep_cache->{translateable_seq})
          if $self->{total_length} && $vep_cache->{translateable_seq};

        $hash->{Protein_position}  = format_coords($tv->translation_start, $tv->translation_end);
        $hash->{Protein_position} .= '/'.length($vep_cache->{peptide})
          if $self->{total_length} && $vep_cache->{peptide};

        $self->add_sift_polyphen($vfoa, $hash);
      }
    }

    # HGVS
    if($self->{hgvs}) {
      my $hgvs_t = $vfoa->hgvs_transcript;
      my $hgvs_p = $vfoa->hgvs_protein;
      my $offset = $vfoa->hgvs_offset;

      # URI encode "="
      $hgvs_p =~ s/\=/\%3D/g if $hgvs_p && !$self->{no_escape};

      $hash->{HGVSc} = $hgvs_t if $hgvs_t;
      $hash->{HGVSp} = $hgvs_p if $hgvs_p;

      $hash->{HGVS_OFFSET} = $offset if $offset;
    }
  }

  return $hash;
}

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

# process RegulatoryFeatureVariationAllele
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

# process MotifFeatureVariationAllele
sub MotifFeatureVariationAllele_to_output_hash {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  # run "super" method
  $hash = $self->VariationFeatureOverlapAllele_to_output_hash(@_);
  return undef unless $hash;

  my $mf = $vfoa->motif_feature;

  # check that the motif has a binding matrix, if not there's not
  # much we can do so don't return anything
  return undef unless defined $mf->binding_matrix;

  my $matrix = ($mf->binding_matrix->description ? $mf->binding_matrix->description.' ' : '').$mf->display_label;
  $matrix =~ s/\s+/\_/g;

  $hash->{Feature_type} = 'MotifFeature';
  $hash->{Feature}      = $mf->binding_matrix->name;
  $hash->{MOTIF_NAME}   = $matrix;
  $hash->{STRAND}       = $mf->strand + 0;
  $hash->{CELL_TYPE}    = $self->get_cell_types($mf) if $self->{cell_type};
  $hash->{MOTIF_POS}    = $vfoa->motif_start if defined $vfoa->motif_start;
  $hash->{HIGH_INF_POS} = ($vfoa->in_informative_position ? 'Y' : 'N');

  my $delta = $vfoa->motif_score_delta if $vfoa->variation_feature_seq =~ /^[ACGT]+$/;
  $hash->{MOTIF_SCORE_CHANGE} = sprintf("%.3f", $delta) if defined $delta;

  return $hash;
}

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

sub IntergenicVariationAllele_to_output_hash {
  my $self = shift;
  return $self->VariationFeatureOverlapAllele_to_output_hash(@_);
}

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
  $hash->{Feature}      = $feature_type eq 'MotifFeature' ? $feature->binding_matrix->name : $feature->stable_id;

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

sub IntergenicStructuralVariationAllele_to_output_hash {
  my $self = shift;
  return $self->BaseStructuralVariationOverlapAllele_to_output_hash(@_);
}

sub plugins {
  return $_[0]->{plugins} || [];
}

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

      $vf->{allele_string} = $vf->{original_allele_string};
      $vf->{start}         = $vf->{original_start};
      $vf->{end}           = $vf->{original_end};

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

        $iv->base_variation_feature($original);

        if(my $oiv = $original->{intergenic_variation}) {
            push @{$oiv->{alt_alleles}}, @{$iv->{alt_alleles}};
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

1;
