=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation - base class for variation (overlap) sources

DO NOT USE DIRECTLY

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  foreach my $existing_vf(@{$self->get_all_features_by_InputBuffer($buffer)}) {
    foreach my $vf(
      grep {ref($_) ne 'Bio::EnsEMBL::Variation::StructuralVariationFeature'}
      @{$buffer->get_overlapping_vfs($existing_vf->{start}, $existing_vf->{end})}
    ) {
      push @{$vf->{existing}}, $existing_vf unless $self->is_var_novel($existing_vf, $vf);
    }
  }

  $self->frequency_check_buffer($buffer) if $self->{check_frequency};
}

sub is_var_novel {
  my $self = shift;
  my $existing_var = shift;
  my $new_var = shift;

  my $is_novel = 1;
  
  my $matched_coords = $existing_var->{start} == $new_var->start && $existing_var->{end} == $new_var->end;
  $is_novel = 0 if $matched_coords;

  # can't compare alleles with e.g. HGMD_MUTATION so just include it
  return 0 if $matched_coords && $existing_var->{allele_string} !~ /\//;

  unless($self->{no_check_alleles}) {
    my %existing_alleles;

    $existing_alleles{$_} = 1 for split '\/', $existing_var->{allele_string};

    my $seen_new = 0;
    foreach my $a(grep {$_ ne 'N'} split '\/', ($new_var->allele_string || "")) {
      reverse_comp(\$a) if $new_var->strand ne $existing_var->{strand};
      $seen_new = 1 unless defined $existing_alleles{$a};
    }

    $is_novel = 1 if $seen_new;
  }

  return $is_novel;
}

sub filter_variation {
  my ($self, $var) = @_;
  return 0 unless $var->{failed} <= (defined($self->{failed}) ? $self->{failed} : 0);
}

sub up_down_size {
  return 1;
}

1;
