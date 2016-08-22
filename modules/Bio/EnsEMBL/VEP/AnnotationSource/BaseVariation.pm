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

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(no_check_alleles failed check_frequency freq_pop freq_freq freq_gt_lt freq_filter)]);

  # add this flag to tell VEP runner to use this first
  $self->{can_filter_vfs} = 1;

  return $self;
}

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

sub check_frequency_filter {
  my $self = shift;

  throw("ERROR: No population specified") unless $self->{freq_pop};

  my $freq_pop_full = uc($self->{freq_pop});
  $freq_pop_full =~ s/EXAC/ExAC/;
  my ($freq_group, $freq_pop) = split('_', $freq_pop_full);

  my %valid = map {$_ => 1} @{$self->get_valid_populations};

  throw("ERROR: Invalid population ".$self->{freq_pop}." specified") unless $freq_pop_full eq '1KG_ALL' || $valid{$freq_pop_full} || $valid{$freq_pop || ''};
}

sub frequency_check_buffer {
  my $self = shift;
  my $buffer = shift;

  my @passed;
  my $freq_filter = $self->{freq_filter};

  foreach my $vf(@{$buffer->buffer}) {
    $self->get_frequency_data($vf);

    # include only
    if($freq_filter eq 'include') {

      # we can disregard this VF straight away if there's no freq data
      next unless $vf->{_freq_check_pass};

      # all passed?
      if($vf->{_freq_check_all_passed}) {
        push @passed, $vf;
      }
      else {
        # new alt alleles
        my @alts = grep {$vf->{_freq_check_pass}->{$_}} @{$vf->alt_alleles};

        # new allele_string
        $vf->allele_string(join('/', $vf->ref_allele_string, @alts));
      }
    }

    # exclude ones that don't pass
    else {

      # no frequency data so we have to include it
      unless($vf->{_freq_check_pass}) {
        push @passed, $vf;
        next;
      }

      # all failed?
      if($vf->{_freq_check_all_failed}) {
        push @passed, $vf;
      }
      elsif(!$vf->{_freq_check_all_passed}) {
        # new alt alleles
        my @alts = grep {!$vf->{_freq_check_pass}->{$_}} @{$vf->alt_alleles};

        # new allele_string
        $vf->allele_string(join('/', $vf->ref_allele_string, @alts));

        push @passed, $vf;
      }
    }
  }

  my $filtered_count = scalar @{$buffer->buffer} - scalar @passed;

  # replace/update the buffer
  if($filtered_count) {

    # we need to reset so that interval trees etc don't get confused
    $buffer->reset_buffer();

    $buffer->buffer(\@passed);

    # update stats
    $self->stats->increment_filtered_variants($filtered_count);
  }
}

sub _add_check_freq_data_to_vf {
  my $self = shift;
  my $vf = shift;
  my $freq_data = shift;

  my $freq_freq     = $self->{freq_freq};
  my $freq_gt_lt    = $self->{freq_gt_lt};
  my $freq_pop_full = uc($self->{freq_pop});

  my $alts = $vf->alt_alleles();

  my %pass;
  my %used_freqs;
  my $pass_count = 0;

  foreach my $alt(@$alts) {
    my $pass = 0;

    # set freq to 0 if it hasn't been observed
    # this I think is a fair assumption if the position has been fully typed
    my $f = $freq_data->{$alt} || 0;

    $pass = 1 if $f >= $freq_freq and $freq_gt_lt eq 'gt';
    $pass = 1 if $f <= $freq_freq and $freq_gt_lt eq 'lt';
    $used_freqs{$alt} = $f;

    $pass_count += $pass;
    $pass{$alt} = $pass;
  }

  $vf->{_freq_check_freqs}->{$freq_pop_full} = \%used_freqs;
  $vf->{_freq_check_pass} = \%pass;

  # all passed or failed?
  $vf->{_freq_check_all_passed} = 1 if $pass_count == scalar @$alts;
  $vf->{_freq_check_all_failed} = 1 if $pass_count == 0;
}

1;
