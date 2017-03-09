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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript - local disk transcript annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::Haplo::AnnotationSource::BaseTranscript Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript);

sub new {  
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->param('haplotype_frequencies', $self->dir.'/haplotype_frequencies.txt.gz') if -e $self->dir.'/haplotype_frequencies.txt.gz';

  return $self;
}

sub _tree_coords_filename {
  return $_[0]->dir.'/transcript_coords.txt';
}

sub _tree_file_data {
  return [$_[1]->seq_region_start, $_[1]->seq_region_end];
}

sub _tree_insert_file_line {
  my ($self, $tree, $line) = @_;
  chomp($line);
  $tree->insert(split("\t", $line));
}

1;