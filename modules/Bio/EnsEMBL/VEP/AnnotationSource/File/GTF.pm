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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::GTF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::GTF - GTF annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::GTF;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::Parser::GTFTabix;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF);

my %PARENTS = (
  'exon' => 'transcript',
  'cds' => 'transcript',
  'transcript' => 'gene',
);

my %INCLUDE_FEATURE_TYPES = map {$_ => 1} qw(
  cds
  CDS
  exon
  gene
  transcript
);

sub parser {
  my $self = shift;
  if(!exists($self->{parser})) {
    $self->{parser} = Bio::EnsEMBL::IO::Parser::GTFTabix->open($self->file);
  }
  return $self->{parser};
  # return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::GTFTabix->open($self->file);
}

sub include_feature_types {
  return \%INCLUDE_FEATURE_TYPES;
}

sub _record_get_parent_id {
  my ($self, $record) = @_;
  if(!exists($record->{_parent_id})) {
    my $type = lc($record->{type});
    $record->{_parent_id} = $PARENTS{$type} ? $record->{attributes}->{$PARENTS{$type}.'_id'} : undef;
  }
  return $record->{_parent_id};
}

sub _record_get_id {
  my ($self, $record) = @_;
  return $record->{_id} ||= $record->{attributes}->{$record->{type}.'_id'} || $record->{md5};
}

sub _record_get_biotype {
  return $_[1]->{attributes}->{transcript_biotype} || $_[1]->{attributes}->{biotype} || $_[1]->{source};
}

1;