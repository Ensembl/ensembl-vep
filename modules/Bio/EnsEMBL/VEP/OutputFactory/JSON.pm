=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
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

# EnsEMBL module for Bio::EnsEMBL::VEP::OutputFactory::JSON
#
#

=head1 NAME

Bio::EnsEMBL::VEP::OutputFactory::JSON - JSON format output factory

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::OutputFactory::JSON;

use base qw(Bio::EnsEMBL::VEP::OutputFactory);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::Constants;

use Bio::EnsEMBL::VEP::Utils qw(numberify);

use JSON;
use Scalar::Util qw(looks_like_number);

my %SKIP_KEYS = (
  'Uploaded_variation' => 1,
  'Location' => 1,
);

my %RENAME_KEYS = (
  'consequence' => 'consequence_terms',
  'gene' => 'gene_id',
  'allele' => 'variant_allele',
  'symbol' => 'gene_symbol',
  'symbol_source' => 'gene_symbol_source',
  'overlapbp' => 'bp_overlap',
  'overlappc' => 'percentage_overlap',
  'refseq' => 'refseq_transcript_ids',
  'ensp' => 'protein_id',
  'chr' => 'seq_region_name',
  'variation_name' => 'id',
);

my %NUMBERIFY_EXEMPT = (
  'seq_region_name' => 1,
  'id' => 1,
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    assembly
    cache_assembly
  )]);

  return $self;
}

sub get_all_lines_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my @return;

  $self->{json_obj} ||= JSON->new;

  foreach my $vf(@{$buffer->buffer}) {

    my $hash = {
      id              => $vf->{variation_name},
      seq_region_name => $vf->{chr},
      start           => $vf->{start},
      end             => $vf->{end},
      strand          => $vf->{strand},
      allele_string   => $vf->{allele_string} || $vf->{class_SO_term},
      assembly_name   => $self->{assembly} || $self->{cache_assembly},
      # _order        => $vf->{_order},
    };

    # add original input for use by POST endpoints
    $hash->{input} = join("\t", @{$vf->{_line}}) if defined($vf->{_line});

    # get other data from super method
    my $extra_hash = $self->VariationFeature_to_output_hash($vf);
    $hash->{lc($_)} = $extra_hash->{$_} for grep {!$SKIP_KEYS{$_}} keys %$extra_hash;

    $self->add_VariationFeatureOverlapAllele_info($vf, $hash);

    $self->add_colocated_variant_info($vf, $hash);

    numberify($hash, \%NUMBERIFY_EXEMPT);
    
    push @return, $self->{json_obj}->encode($hash);
  }

  return \@return;
}

sub add_VariationFeatureOverlapAllele_info {
  my $self = shift;
  my $vf = shift;
  my $hash = shift;

  # record all cons terms so we can get the most severe
  my @con_terms;

  # add consequence stuff
  foreach my $vfoa_hash(@{$self->get_all_VariationFeatureOverlapAllele_output_hashes($vf, {})}) {

    # lc and remove empty
    foreach my $key(keys %$vfoa_hash) {
      my $tmp = $vfoa_hash->{$key};
      delete $vfoa_hash->{$key};

      next if !defined($tmp) || $tmp eq '-';

      # convert YES to 1
      $tmp = 1 if $tmp eq 'YES';

      # fix position fields into start and end
      if($key =~ /(\w+?)\_position$/i) {
        my $coord_type = lc($1);
        my ($s, $e) = split('-', $tmp);
        $vfoa_hash->{$coord_type.'_start'} = $s;
        $vfoa_hash->{$coord_type.'_end'} = defined($e) && $e =~ /^\d+$/ ? $e : $s;

        # on rare occasions coord can be "?"; for now just don't print anything
        delete $vfoa_hash->{$coord_type.'_start'} unless looks_like_number($vfoa_hash->{$coord_type.'_start'});
        delete $vfoa_hash->{$coord_type.'_end'}   unless looks_like_number($vfoa_hash->{$coord_type.'_end'});
        next;
      }

      $vfoa_hash->{lc($key)} = $tmp;
    }

    my $ftype = lc($vfoa_hash->{feature_type} || 'intergenic');
    $ftype =~ s/feature/\_feature/;
    delete $vfoa_hash->{feature_type};

    # fix SIFT and PolyPhen
    foreach my $tool(qw(sift polyphen)) {
      if(defined($vfoa_hash->{$tool}) && $vfoa_hash->{$tool} =~ m/([a-z\_]+)?\(?([\d\.]+)?\)?/i) {
        my ($pred, $score) = ($1, $2);
        $vfoa_hash->{$tool.'_prediction'} = $pred if $pred;
        $vfoa_hash->{$tool.'_score'} = $score if defined($score);
        delete $vfoa_hash->{$tool};
      }
    }

    # fix domains
    if(defined($vfoa_hash->{domains})) {
      my @dom;

      foreach(@{$vfoa_hash->{domains}}) {
        m/(\w+)\:(\w+)/;
        push @dom, {"db" => $1, "name" => $2} if $1 && $2;
      }
      $vfoa_hash->{domains} = \@dom;
    }

    # log observed consequence terms
    push @con_terms, @{$vfoa_hash->{consequence}};

    # rename
    my %rename = %RENAME_KEYS;

    $rename{feature} = lc($ftype).'_id';
    foreach my $key(grep {defined($vfoa_hash->{$_})} keys %rename) {
      $vfoa_hash->{$rename{$key}} = $vfoa_hash->{$key};
      delete $vfoa_hash->{$key};
    }

    push @{$hash->{$ftype.'_consequences'}}, $vfoa_hash;
  }

  my %all_cons = %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
  $hash->{most_severe_consequence} = (sort {$all_cons{$a}->rank <=> $all_cons{$b}->rank} grep {$_ ne '?'} @con_terms)[0] || '?';

  return $hash;
}

sub add_colocated_variant_info {
  my $self = shift;
  my $vf = shift;
  my $hash = shift;

  foreach my $ex_orig(@{$vf->{existing} || []}) {

    # work on a copy as we're going to modify/delete things
    my $ex;
    %$ex = %$ex_orig;

    delete $ex->{$_} for qw(failed);

    # frequencies
    foreach my $pop(grep {defined($ex->{$_})} qw(
      AFR AMR ASN EAS SAS EUR
      AA EA
      ExAC ExAC_Adj ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS
    )) {
      my $tmp = $ex->{$pop};
      my $lc_pop = lc($pop);

      if($tmp =~ /(\w)\:([\d\.\-e]+)/) {
        $ex->{$lc_pop.'_maf'} = $2;
        $ex->{$lc_pop.'_allele'} = $1;
      }
      else {
        $ex->{$lc_pop.'_maf'} = $tmp;
      }

      delete $ex->{$pop};
    }

    # remove empty
    foreach my $key(keys %$ex) {
      delete $ex->{$key} if !defined($ex->{$key}) || $ex->{$key} eq '' || ($key !~ /maf/ && $ex->{$key} eq 0);
    }

    # rename
    foreach my $key(grep {defined($ex->{$_})} keys %RENAME_KEYS) {
      $ex->{$RENAME_KEYS{$key}} = $ex->{$key};
      delete $ex->{$key};
    }

    push @{$hash->{colocated_variants}}, $ex;
  }

  return $hash;
}

1;
