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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat - local disk regulatory feature annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use parent qw(Bio::EnsEMBL::VEP::AnnotationSource::Cache);

sub get_dump_file_name {
  my $self = shift;
  my $chr  = shift;
  my $region = shift;

  throw("No chromosome given") unless $chr;
  throw("No region given") unless $region;

  # allow to pass region (start-end) or $start, $end
  $region .= '-'.shift if @_;

  return sprintf(
    "%s/%s/%s_reg\.%s",
    $self->root_dir,
    $chr,
    $region,
    $self->file_suffix
  );
}

1;