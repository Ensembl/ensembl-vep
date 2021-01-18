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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache - local disk annotation source

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Base class for all cache-based AnnotationSource classes.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);


=head2 dir

  Example    : $dir = $as->dir();
  Description: Gets the path to the directory containing the cache data
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub dir {
  my $self = shift;
  $self->{dir} = shift if @_;
  return $self->{dir};
}


=head2 valid_chromosomes

  Example    : $chrs = $as->valid_chromosomes();
  Description: Gets valid chromosome names for this cache
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub valid_chromosomes {
  return $_[0]->{valid_chromosomes} || [];
}

1;