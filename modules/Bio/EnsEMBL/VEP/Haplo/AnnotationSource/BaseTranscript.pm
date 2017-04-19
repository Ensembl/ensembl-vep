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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::AnnotationSource::BaseTranscript - base haplotype annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationSource::BaseTranscript;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::SampleGenotypeFeature;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my @return;

  my $samples = $buffer->parser->samples;

  foreach my $tr(@{$self->get_all_features_by_InputBuffer($buffer)}) {
    my ($s, $e);
    if(ref($tr) eq 'HASH') {
      ($s, $e) = ($tr->{start}, $tr->{end});
    }
    else {
      ($s, $e) = ($tr->seq_region_start, $tr->seq_region_end);
    }

    my $vfs = $buffer->get_overlapping_vfs($s, $e);
    next unless @$vfs;

    $tr = $self->lazy_load_transcript($tr);
    next unless $tr && $tr->{biotype} eq 'protein_coding';

    next if $self->filter_set && !$self->filter_set->evaluate($tr);

    my @gts;

    foreach my $vf_hash(@$vfs) {
      push @gts, @{$self->get_genotypes($vf_hash, $samples)};
    }

    if(@gts) {
      my $ct = $self->create_container($tr, \@gts, [values %$samples]);
      push @return, $ct if $ct;
    }
  }

  return \@return;
}

sub create_vf {
  my ($self, $hash) = @_;

  my @alleles = split(',', $hash->{alleles});

  return Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => $hash->{start},
    end            => $hash->{end},
    allele_string  => join('/', @alleles),
    strand         => 1,
    map_weight     => 1,
    adaptor        => $self->get_adaptor('variation', 'VariationFeature'),
    variation_name => $hash->{ids}->[0] || undef,
    chr            => $hash->{chr},
    slice          => $self->get_slice($hash->{chr}),
  });
}

sub get_genotypes {
  my ($self, $hash, $samples) = @_;

  if(!exists($hash->{gt_objects})) {
    my @gts;
    my $vf = $hash->{vf} ||= $self->create_vf($hash);

    foreach my $sample(keys %{$hash->{gts}}) {
      push @gts, Bio::EnsEMBL::Variation::SampleGenotypeFeature->new_fast({
        # _variation_id     => $vf->{_variation_id},
        variation_feature => $vf,
        sample            => $samples->{$sample},
        genotype          => [split(/\||\/|\\/, $hash->{gts}->{$sample})],
        phased            => 1,
        start             => $vf->start,
        end               => $vf->end,
        strand            => $vf->seq_region_strand,
        slice             => $vf->slice,
      })
    }

    $hash->{gt_objects} = \@gts;
  }

  return $hash->{gt_objects};
}

sub create_container {
  my ($self, $tr, $gts, $samples) = @_;

  return Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
    -transcript => $tr,
    -genotypes  => $gts,
    -samples    => $samples,
  );
}

1;