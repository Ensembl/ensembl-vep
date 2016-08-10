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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::GFF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::GFF - GFF annotation source

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
  lincRNA
  lincRNA_gene
  miRNA
  miRNA_gene
  mt_gene
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

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::GFF3Tabix->open($self->file);
}

sub include_feature_types {
  return \%INCLUDE_FEATURE_TYPES;
}

1;