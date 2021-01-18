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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::GTF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::GTF - GTF annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GTF->new({
  config => $config,
  file   => "my_genes.gtf.gz"
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

GTF format custom annotation source. GTF files must be chromosome/pos
sorted, compressed with bgzip and indexed with tabix.

Other aspects/limitations discussed here:

http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff

=head1 METHODS

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
  'stop_codon' => 'transcript',
  'transcript' => 'gene',
);

my %INCLUDE_FEATURE_TYPES = map {$_ => 1} qw(
  cds
  CDS
  stop_codon
  exon
  gene
  transcript
);


=head2 parser

  Example    : $parser = $as->parser();
  Description: Get ensembl-io parser to read from file
  Returntype : Bio::EnsEMBL::IO::Parser::GTFTabix
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::GTFTabix->open($self->file, must_parse_metadata => 0);
}


=head2 include_feature_types

  Example    : $types = $as->include_feature_types();
  Description: Get hashref of GFF record types this class can parse
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub include_feature_types {
  return \%INCLUDE_FEATURE_TYPES;
}


=head2 _record_get_parent_ids

  Arg 1      : hashref $record_hash
  Example    : $ids = $as->_record_get_parent_ids($record);
  Description: Get IDs of parent records for this record
  Returntype : listref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _record_get_parent_ids {
  my ($self, $record) = @_;
  if(!exists($record->{_parent_id})) {
    my $type = lc($record->{type});
    my @parent_ids = ($PARENTS{$type}) ? split(',',$record->{attributes}->{$PARENTS{$type}.'_id'}) : ();
    $record->{_parent_id} = \@parent_ids;
  }
  return $record->{_parent_id};
}


=head2 _record_get_id

  Arg 1      : hashref $record_hash
  Example    : $id = $as->_record_get_id($record);
  Description: Get ID of this record
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _record_get_id {
  my ($self, $record) = @_;
  return $record->{_id} ||= $record->{attributes}->{$record->{type}.'_id'} || $record->{md5};
}


=head2 _record_get_biotype

  Arg 1      : hashref $transcript_record_hash
  Example    : $biotype = $as->_record_get_biotype($tr_record);
  Description: Get sequence ontology (SO) biotype of this record. Attempts to
               find it in the "biotype", "transcript_type" or "transcript_biotype"
               attribute fields, and if that fails default to the source field.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _record_get_biotype {
  return
    $_[1]->{attributes}->{transcript_biotype} ||
    $_[1]->{attributes}->{transcript_type} ||
    $_[1]->{attributes}->{biotype} ||
    $_[1]->{source};
}

1;
