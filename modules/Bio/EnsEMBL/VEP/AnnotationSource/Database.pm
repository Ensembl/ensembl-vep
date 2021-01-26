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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database - database annotation source

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Base class for all database-based AnnotationSource classes.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config $config,
                 filter => string $transcript_filter,
                 bam    => string $bam_file,
               }
  Example    : Not invoked directly
  Description: Creates a new database AnnotationSource object
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::Database
  Exceptions : none
  Caller     : AnnotationSourceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->{cache_region_size} ||= 50000;

  return $self;
}


=head2 get_slice

  Arg 1      : string $chr
  Example    : $slice = $obj->get_slice('MT')
  Description: Gets whole-chromosome slices from slice adaptor.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub get_slice {
  my $self = shift;
  my $chr = shift;
  my $chr_is_seq_region = shift;

  if(!exists($self->{_slice_cache}->{$chr})) {
    my $sa = $self->get_adaptor(
      $self->{core_type} || $self->param('core_type'),
      'Slice'
    );

    my $slice = $chr_is_seq_region ? $sa->fetch_by_seq_region_id($chr) : $sa->fetch_by_region(undef, $chr);
    $slice->is_circular if $slice;

    $self->{_slice_cache}->{$chr} = $slice;
  }

  return $self->{_slice_cache}->{$chr};
}


=head2 valid_chromosomes

  Example    : $chrs = $as->valid_chromosomes();
  Description: Gets valid chromosome names for this database
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub valid_chromosomes {
  return $_[0]->{valid_chromosomes} ||= [sort keys %{$_[0]->chr_lengths}];
}

1;