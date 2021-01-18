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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Database::Transcript - database transcript annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
  config => $config
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

AnnotationSource that reads transcript data from a core schema database.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Database::Transcript;

use base qw(Bio::EnsEMBL::VEP::Haplo::AnnotationType::Transcript Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript);


=head2 populate_tree

  Arg 1      : Bio::EnsEMBL::VEP::TranscriptTree
  Example    : $as->populate_tree($tree);
  Description: Populates the given transcript tree with data from the transcripts
               in the database.
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub populate_tree {
  my ($self, $tree) = @_;
  my $ta = $self->get_adaptor('core', 'transcript');
  $tree->insert($_->seq_region_name, $_->seq_region_start, $_->seq_region_end) for @{$ta->fetch_all_by_biotype('protein_coding')};
}


=head2 create_container

  Arg 1      : Bio::EnsEMBL::Transcript
  Arg 2      : arrayref of Bio::EnsEMBL::Variation::SampleGenotypeFeature
  Arg 3      : arrayref of Bio::EnsEMBL::Variation::Sample
  Example    : $thc = $as->create_container($tr, $gts, $samples);
  Description: Creates a TranscriptHaplotypeContainer for the given transcript
  Returntype : Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::Haplo::BaseTranscript
  Status     : Stable

=cut

sub create_container {
  my ($self, $tr, $gts, $samples) = @_;

  return Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
    -transcript => $tr,
    -genotypes  => $gts,
    -samples    => $samples,
    -db         => $self->get_adaptor('variation', 'variation')->db,
  );
}

1;