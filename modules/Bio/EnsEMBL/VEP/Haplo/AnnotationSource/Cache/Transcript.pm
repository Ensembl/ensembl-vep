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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript - local disk transcript annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
  config => $config,
  dir => $dir,
  cache_region_size => 1e6,
  info => $version_data,
  valid_chromosomes => $valid_chromosomes,
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

AnnotationSource that reads transcript data from a VEP cache.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::Haplo::AnnotationType::Transcript Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript);


=head2 new

  Arg 1      : hashref $args
               {
                 config            => Bio::EnsEMBL::VEP::Config,
                 dir               => string,
                 cache_region_size => int (1e6),
                 info              => hashref,
                 valid_chromosomes => listref of strings,
                 serializer_type   => string (storable, sereal),
                 filter            => string,
                 bam               => string,
               }
  Example    : $ib = Bio::EnsEMBL::VEP::Haplo::InputBuffer->new($args);
  Description: Creates a Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript
               object. Sets the haplotype_frequencies param if an appropriate file
               is found in dir.
  Returntype : Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::CacheDir
  Status     : Stable

=cut

sub new {  
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->param('haplotype_frequencies', $self->dir.'/haplotype_frequencies.txt.gz') if -e $self->dir.'/haplotype_frequencies.txt.gz';

  return $self;
}


=head2 _tree_coords_filename

  Example    : $file = $as->_tree_coords_filename();
  Description: Gets the filename for the transcript tree
  Returntype : string
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _tree_coords_filename {
  return $_[0]->dir.'/transcript_coords.txt';
}


=head2 _tree_file_data

  Arg 1      : Bio::EnsEMBL::Transcript $transcript
  Example    : $data = $as->_tree_file_data($transcript);
  Description: Gets the data required to write to the transcript tree file
               from a given transcript.
  Returntype : arrayref [$start, $end]
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _tree_file_data {
  return [$_[1]->seq_region_start, $_[1]->seq_region_end];
}


=head2 _tree_insert_file_line

  Arg 1      : Bio::EnsEMBL::VEP::TranscriptTree $tree
  Arg 2      : string $tab_delimited_data
  Example    : $as->_tree_insert_file_line($tree, $data);
  Description: Inserts into a transcript tree a line of data read from
               a transcript tree file
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _tree_insert_file_line {
  my ($self, $tree, $line) = @_;
  chomp($line);
  $tree->insert(split("\t", $line));
}

1;