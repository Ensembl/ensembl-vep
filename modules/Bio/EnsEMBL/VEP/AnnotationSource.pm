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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource - Base class used for VEP annotation sources

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

sub get_all_features_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my $regions = $self->get_all_regions_by_InputBuffer($buffer);

  my @not_cached;
  my @features;

  # attempt to fetch from memory first
  foreach my $region(@$regions) {
    my $cache_features = $self->get_features_by_regions_cached([$region]);

    # if successful, add them to the @features to be returned
    if(scalar @$cache_features) {
      push @features, @$cache_features;
    }

    # otherwise keep track of this region as being unloaded
    else {
      push @not_cached, $region;
    }
  }

  # now get the remaining unloaded ones from disk
  push @features, @{$self->get_features_by_regions_uncached(\@not_cached)};

  $self->clean_cache($regions);

  # use merge_features to remove duplicates
  return $self->merge_features(\@features);
}

sub get_all_regions_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;
  my $cache_region_size = shift || $self->{cache_region_size};

  my @regions = ();
  my %seen = ();
  my $up_down_size = $self->up_down_size();

  foreach my $vf(@{$buffer->buffer}) {
    my $chr = $vf->{chr} || $vf->slice->seq_region_name;
    throw("ERROR: Cannot get chromosome from VariationFeature") unless $chr;

    my @region_starts =
      $vf->{start} > $vf->{end} ?
      ($vf->{end}   - $up_down_size, $vf->{start} + $up_down_size) :
      ($vf->{start} - $up_down_size, $vf->{end}   + $up_down_size);

    foreach my $region_start(map {int($_ / $cache_region_size)} @region_starts) {
      my $key = join(':', ($chr, $region_start));
      next if $seen{$key};

      push @regions, [$chr, $region_start];
      $seen{$key} = 1;
    }
  }

  return \@regions;
}

sub get_features_by_regions_cached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;

  return [map {@{$cache->{$_->[0]}->{$_->[1]} || []}} @$regions];
}

sub cache {
  return $_[0]->{_cache} ||= {};
}

# used after each disk fetch
# removes out-of-range features from cache
sub clean_cache {
  my $self = shift;
  my $keep_regions = shift || [];

  my $cache = $self->cache;

  # copy the regions to be kept to has structure equivalent to $cache
  my $keep_regions_hash;
  $keep_regions_hash->{$_->[0]}->{$_->[1]} = 1 for @$keep_regions;

  # now go through and delete what we don't need any more
  foreach my $chr(keys %$cache) {
    if($keep_regions_hash->{$chr}) {
      delete $cache->{$chr}->{$_} for grep {!$keep_regions_hash->{$chr}->{$_}} keys %{$cache->{$chr}};
    }
    else {
      delete $cache->{$chr};
    }
  }
}

1;