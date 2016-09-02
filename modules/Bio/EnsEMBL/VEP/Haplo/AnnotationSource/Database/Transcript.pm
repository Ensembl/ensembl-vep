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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Database::Transcript - database transcript annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Database::Transcript;

use base qw(Bio::EnsEMBL::VEP::Haplo::AnnotationSource::BaseTranscript Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript);

sub populate_tree {
  my ($self, $tree) = @_;
  my $ta = $self->get_adaptor('core', 'transcript');
  $tree->insert($_->seq_region_name, $_->seq_region_start, $_->seq_region_end) for @{$ta->fetch_all_by_biotype('protein_coding')};
}

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