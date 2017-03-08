=head1 LICENSE

Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->{cache_region_size} ||= 50000;

  return $self;
}

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

sub valid_chromosomes {
  my $self = shift;

  if(!exists($self->{valid_chromosomes})) {
    my $sa = $self->get_adaptor(
      $self->{core_type} || $self->param('core_type'),
      'Slice'
    );

    my @valid_chromosomes;
    my %chr_lengths;

    foreach my $slice(@{$sa->fetch_all('toplevel')}) {
      my $chr = $slice->seq_region_name;
      push @valid_chromosomes, $chr;
      $chr_lengths{$chr} = $slice->length;
    }

    $self->{valid_chromosomes} = \@valid_chromosomes;
    $self->stats->log_db_chromosomes(\%chr_lengths);
  }

  return $self->{valid_chromosomes};
}

1;