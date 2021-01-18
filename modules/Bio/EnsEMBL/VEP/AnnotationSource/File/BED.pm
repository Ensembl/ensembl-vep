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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::BED
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::BED - BED annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::BED->new({
  config => $config,
  file   => "my_features.bed.gz",
  type   => "overlap"
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

BED format custom annotation source. BED files must be chromosome/pos
sorted, compressed with bgzip and indexed with tabix.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::BED;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::Parser::BedTabix;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File);


=head2 parser

  Example    : $parser = $as->parser();
  Description: Get ensembl-io parser to read from file
  Returntype : Bio::EnsEMBL::IO::Parser::BedTabix
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::BedTabix->open($self->file);
}

1;