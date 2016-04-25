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

use Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);

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
}

# gets variation cache columns
sub get_cache_columns {
  my $self = shift;

  if(!exists($self->{cols})) {
    my @copy = @Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation::VAR_CACHE_COLS;
    $self->{cols} = \@copy;
    # push @{$self->{cols}}, 'pubmed' if $self->have_pubmed() && $self->param('pubmed');
    # push @{$self->{cols}}, @{$self->{freq_file_pops}} if defined($self->{freq_file_pops});
  }

  return $self->{cols};
}

sub is_var_novel {
  my $self = shift;
  my $existing_var = shift;
  my $new_var = shift;

  my $is_novel = 1;

  $is_novel = 0 if $existing_var->{start} == $new_var->start && $existing_var->{end} == $new_var->end;

  if($self->{check_alleles}) {
    my %existing_alleles;

    $existing_alleles{$_} = 1 for split '\/', $existing_var->{allele_string};

    my $seen_new = 0;
    foreach my $a(split '\/', ($new_var->allele_string || "")) {
      reverse_comp(\$a) if $new_var->strand ne $existing_var->{strand};
      $seen_new = 1 unless defined $existing_alleles{$a};
    }

    $is_novel = 1 if $seen_new;
  }

  return $is_novel;
}

sub filter_variation {
  my ($self, $var) = @_;
  return 0 unless $var->{failed} <= ($self->{failed} // 0);
}

sub up_down_size {
  return 1;
}

1;
