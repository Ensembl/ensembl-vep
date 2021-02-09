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

Bio::EnsEMBL::VEP::Haplo::AnnotationSource::File::GTF - GTF transcript annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GTF->new({
  config     => $config,
  file       => $gtf_file,
  short_name => $short_name,
  type       => 'overlap',
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

AnnotationSource that reads transcript data from a tabix-indexed GTF file.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationSource::File::GTF;

use base qw(Bio::EnsEMBL::VEP::Haplo::AnnotationType::Transcript Bio::EnsEMBL::VEP::AnnotationSource::File::GTF);


=head2 populate_tree

  Arg 1      : Bio::EnsEMBL::VEP::TranscriptTree
  Example    : $as->populate_tree($tree);
  Description: Populates the given transcript tree with data from the transcripts
               in the GTF file.
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub populate_tree {
  my ($self, $tree) = @_;

  my $parser = $self->parser;

  my %include = %{$self->include_feature_types};
  delete $include{$_} for qw(exon CDS stop_codon);

  foreach my $chr(@{$parser->{tabix_file}->seqnames}) {
    $parser->seek($chr, 1, 1e10);
    while($parser->next) {
      $tree->insert($chr, $parser->get_start, $parser->get_end) if $include{$parser->get_type};
    }
  }

  delete $self->{parser};
}

1;