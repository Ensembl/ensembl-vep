=head1 LICENSE

Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationType::Variation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationType::Variation - base class for variation (overlap) sources

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Helper class for all variant-based AnnotationSource classes. 

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationType::Variation;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $as->annotate_InputBuffer($ib);
  Description: Gets overlapping known variants for the variants in
               the input buffer, checks if they are "novel" (see is_var_novel())
               and adds them to the key "existing" on the VariationFeature
               object. May also filter the input buffer by comparing frequencies
               of known variants to user-specified thresholds.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

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


=head2 get_cache_columns

  Example    : my $cols = $as->get_cache_columns();
  Description: Gets the default columns for the known variant cache file.
               In newer caches the columns are specified in the cache's
               info.txt file, so this is generally only used when creating
               caches.
  Returntype : none
  Exceptions : none
  Caller     : parse_variation(), DumpVEP pipeline
  Status     : Stable

=cut

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


=head2 is_var_novel

  Arg 1      : hashref $known_var
  Arg 2      : Bio::EnsEMBL::Variation::VariationFeature $input_vf
  Example    : $is_novel = $as->is_var_novel($known_var, $input_vf);
  Description: Compares a VariationFeature as created from user input to
               a known variant from the cache/database and determines if
               the input variant is novel. To do this the alleles of each
               are compared - if any of the alleles in the input variant
               are not found in the known variant (accounting for strand)
               then the input variant is considered to be novel.
  Returntype : bool
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

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


=head2 filter_variation

  Arg 1      : hashref $known_var
  Example    : $include = $as->filter_variation($known_var);
  Description: Filter known variant as loaded from the cache/database
               according to its failed (QC) state by comparing to the
               VEP default or user-set state.
  Returntype : bool
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub filter_variation {
  my ($self, $var) = @_;
  return 0 unless $var->{failed} <= (defined($self->{failed}) ? $self->{failed} : 0);
}


=head2 up_down_size

  Example    : $size = $as->up_down_size();
  Description: Gets range in bp that should be added to boundaries
               when fetching features. 1 for variants to allow for
               indel oddness.
  Returntype : int
  Exceptions : none
  Caller     : get_all_regions_by_InputBuffer()
  Status     : Stable

=cut

sub up_down_size {
  return 1;
}

1;
