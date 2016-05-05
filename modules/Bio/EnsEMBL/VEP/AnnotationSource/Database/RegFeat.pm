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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat - database RegFeat annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Database
  Bio::EnsEMBL::VEP::AnnotationSource::BaseRegFeat
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    cell_type
  )]);

  return $self;
}

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;
  my @return;

  my $cache_region_size = $self->{cache_region_size};

  foreach my $region(@{$regions}) {
    my ($c, $region_start) = @$region;

    my $slice = $self->get_slice($c);

    next unless $slice;

    # get a seq_region_Slice as for patch regions $slice won't cover the whole seq_region
    my $sr_slice = $slice->seq_region_Slice();

    my ($s, $e) = map {($_ - $slice->start) + 1} (
      ($region_start * $cache_region_size) + 1,
      ($region_start + 1) * $cache_region_size
    );

    # sanity check start and end
    $s = 1 if $s < 1;
    $e = $slice->length if $e > $slice->length;

    # get sub-slice
    my $sub_slice = $slice->sub_Slice($s, $e);

    next unless $sub_slice;

    my @region_features;

    foreach my $type(qw(RegulatoryFeature MotifFeature)) {
      my $features = $self->get_adaptor('funcgen', $type)->fetch_all_by_Slice($sub_slice);
      next unless defined($features);

      # cell types
      if($self->{cell_type} && scalar(@{$self->{cell_type}})) {
        foreach my $rf(@$features) {

          my %cl;

          # get cell type by fetching all from stable ID
          if($type eq 'RegulatoryFeature') {
            %cl = map {
              $_->feature_set->cell_type->name => $_->feature_type->name
            } @{$rf->adaptor->fetch_all_by_stable_ID($rf->stable_id)};
          }

          # get cell type by fetching regfeats that contain this MotifFeature
          elsif($type eq 'MotifFeature') {
            %cl = map {
              $_->feature_set->cell_type->name => $_->feature_type->name
            } @{$self->get_adaptor('funcgen', 'RegulatoryFeature')->fetch_all_by_attribute_feature($rf)};
          }

          $rf->{cell_types} = \%cl;
        }
      }

      push @region_features,
        map { $_->{_vep_feature_type} ||= $type; $_ }
        map { $_->transfer($sr_slice) }
        @{$features};
    }

    $cache->{$c}->{$region_start} = \@region_features;

    push @return, @region_features;
  }

  return \@return;
}

1;