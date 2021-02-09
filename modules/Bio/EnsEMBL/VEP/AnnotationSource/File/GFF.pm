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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::GFF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::GFF - GFF annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({
  config => $config,
  file   => "my_genes.gff.gz"
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

GFF3 format custom annotation source. GFF files must be chromosome/pos
sorted, compressed with bgzip and indexed with tabix.

Other aspects/limitations discussed here:

http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::GFF;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::Parser::GFF3Tabix;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF);

my %INCLUDE_FEATURE_TYPES = map {$_ => 1} qw(
  aberrant_processed_transcript
  CDS
  C_gene_segment
  exon
  gene
  J_gene_segment
  lnc_RNA
  lincRNA
  lincRNA_gene
  miRNA
  miRNA_gene
  mt_gene
  nc_primary_transcript
  NMD_transcript_variant
  processed_pseudogene
  processed_transcript
  pseudogene
  pseudogenic_transcript
  RNA
  rRNA
  rRNA_gene
  snoRNA
  snoRNA_gene
  snRNA
  snRNA_gene
  supercontig
  transcript
  VD_gene_segment
  V_gene_segment

  CDS
  C_gene_segment
  D_gene_segment
  exon
  gene
  J_gene_segment
  mRNA
  ncRNA
  primary_transcript
  rRNA
  transcript
  tRNA
  V_gene_segment
);


=head2 parser

  Example    : $parser = $as->parser();
  Description: Get ensembl-io parser to read from file
  Returntype : Bio::EnsEMBL::IO::Parser::GFF3Tabix
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::GFF3Tabix->open($self->file, must_parse_metadata => 0);
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
    my $attributes = $record->{attributes};

    my $parents = $attributes->{Parent} || $attributes->{parent};
    # An record can be linked to more than one parent records
    # e.g. an exon can be linked to several transcripts
    my @parent_ids = ($parents) ? split(',',$parents) : ();
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

  if(!exists($record->{_id})) {
    my $attributes = $record->{attributes};
    $record->{_id} = $attributes->{ID} || $attributes->{Name} || $attributes->{id} || $attributes->{name};
  }
  # Use MD5 as ID if it's missing in the record
  $record->{_id} = $record->{md5} if(!$record->{_id} || $record->{_id} !~ /\w+/);

  return $record->{_id};
}


=head2 _record_get_biotype

  Arg 1      : hashref $transcript_record_hash
  Arg 2      : hashref $gene_record_hash
  Example    : $biotype = $as->_record_get_biotype($tr_record, $gene_record);
  Description: Get sequence ontology (SO) biotype of this record. Attempts to
               find it in the "biotype" or "transcript_type" attribute fields,
               and if that fails (as it will for RefSeq GFFs), make an
               educated guess looking at the record type and various other
               attributes.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub _record_get_biotype {
  my ($self, $record, $gene_record) = @_;

  if(!exists($record->{_biotype})) {

    # Ensembl-y GFFs have biotype as an attribute
    my $biotype = $record->{attributes}->{biotype} || $record->{attributes}->{transcript_type};

    # others we need to (guess) work it out
    if(!$biotype) {
      my $type = lc($record->{type});

      if($type eq 'mrna') {
        $biotype = 'protein_coding';
      }
      elsif($type eq 'ncrna') {
        $biotype = $record->{attributes}->{ncrna_class};
      }
      elsif($type =~ /^([a-z]+)_gene_segment$/) {
        $biotype = 'IG_'.uc($1).'_gene';
      }
      elsif($gene_record && ($gene_record->{attributes}->{description} || '') =~ /^microRNA/) {
        $biotype = 'miRNA';
      }
      elsif($record->{attributes}->{gbkey}) {
        $biotype = $record->{attributes}->{gbkey};
      }
    }

    $record->{_biotype} = $biotype;
  }

  return $record->{_biotype};
}

1;
