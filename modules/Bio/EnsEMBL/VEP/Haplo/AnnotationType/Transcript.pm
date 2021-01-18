=head1 LICENSE

Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::Haplo::AnnotationType::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::AnnotationType::Transcript - base haplotype annotation source

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Helper class for all Haplo::AnnotationSource classes. Contains the bulk of the
code that carries out the annotation.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationType::Transcript;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::SampleGenotypeFeature;
use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::Haplo::InputBuffer
  Example    : my $containers = $as->annotate_InputBuffer($ib);
  Description: Creates TranscriptHaplotypeContainers for the variants in
               the InputBuffer for each Transcript that they overlap.
  Returntype : arrayref of Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::Haplo::Runner
  Status     : Stable

=cut

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my @return;

  my $samples = $buffer->parser->samples;
  my $io_parser = $buffer->parser->parser;

  foreach my $tr(@{$self->get_all_features_by_InputBuffer($buffer)}) {

    $tr = $self->lazy_load_transcript($tr);
    next unless $tr && $tr->{biotype} eq 'protein_coding';

    next if $self->filter_set && !$self->filter_set->evaluate($tr);

    # we only want variants overlapping the exons
    my @vfs = 
      map {@{$buffer->get_overlapping_vfs($_->seq_region_start, $_->seq_region_end)}}
      @{$tr->{_variation_effect_feature_cache}->{sorted_exons} || $tr->get_all_Exons};
    next unless @vfs;

    my @gts;

    foreach my $vf_hash(@vfs) {
      push @gts, @{$self->get_genotypes($vf_hash, $samples, $io_parser)};
    }

    if(@gts) {
      my $ct = $self->create_container($tr, \@gts, [values %$samples]);
      push @return, $ct if $ct;
    }
  }

  return \@return;
}


=head2 create_vf

  Arg 1      : hashref $vf_hash
  Example    : my $vf = $as->create_vf($vf_hash);
  Description: Creates VariationFeature objects from the hashes that are
               created by the InputBuffer's parser.
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Caller     : get_genotypes()
  Status     : Stable

=cut

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


=head2 get_genotypes

  Arg 1      : hashref $vf_hash
  Arg 2      : arrayref of Bio::EnsEMBL::Variation::Sample $samples
  Arg 3      : Bio::EnsEMBL::IO::Parser::VCF4 $parser
  Example    : my $gts = $as->get_genotypes($vf_hash, $samples, $parser);
  Description: Creates SampleGenotypeFeature objects from the variant hashes
               and samples from the InputBuffer. The parser object is used
               to lazy-load genotypes from the VCF record used to generate
               the $vf_hash.
  Returntype : arrayref of Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub get_genotypes {
  my ($self, $hash, $samples, $parser) = @_;

  if(!exists($hash->{gt_objects})) {
    my @gts;
    my $vf = $hash->{vf} ||= $self->create_vf($hash);

    # "re-inject" the raw VCF record into the parser object
    # this allows us to lazy-fetch the genotypes here
    my $bak = $parser->{record};
    $parser->{record} = $hash->{record};
    my $raw_gts = $parser->get_samples_genotypes(undef, 1);
    $parser->{record} = $bak;

    foreach my $sample(keys %$raw_gts) {
      push @gts, Bio::EnsEMBL::Variation::SampleGenotypeFeature->new_fast({
        # _variation_id     => $vf->{_variation_id},
        variation_feature => $vf,
        sample            => $samples->{$sample},
        genotype          => [split(/\||\/|\\/, $raw_gts->{$sample})],
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


=head2 create_container

  Arg 1      : Bio::EnsEMBL::Transcript
  Arg 2      : arrayref of Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Arg 3      : arrayref of Bio::EnsEMBL::Variation::Sample
  Example    : $thc = $as->create_container($tr, $gts, $samples);
  Description: Creates a TranscriptHaplotypeContainer for the given transcript
  Returntype : Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub create_container {
  my ($self, $tr, $gts, $samples) = @_;

  return Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
    -transcript => $tr,
    -genotypes  => $gts,
    -samples    => $samples,
  );
}

1;