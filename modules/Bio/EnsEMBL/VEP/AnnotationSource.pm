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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource - Base class used for VEP annotation sources

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Base class for all VEP annotation sources. Contains basic methods for retrieving
data from the annotation source; some or most of these may be overridden in any
given child class.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::VEP::FilterSet;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);


=head2 new

  Arg 1      : hashref $args
  Example    : Invoked by child classes only
  Description: Creates a new AnnotationSource. Copies all values in $args to $self.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource
  Exceptions : none
  Caller     : Child classes
  Status     : Stable

=cut

sub new {  
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  $self->{$_} = $hashref->{$_} for keys %$hashref;

  if(my $filter = $hashref->{filter}) {
    my @set = ref($filter) eq 'ARRAY' ? @$filter : ($filter);
    $self->filter_set(Bio::EnsEMBL::VEP::FilterSet->new(@set)) if scalar @set;
  }

  return $self;
}


=head2 get_all_features_by_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Example    : $features = $as->get_all_features_by_InputBuffer($ib)
  Description: Fetches all features that overlap the variants currently
               in the input buffer.
  Returntype : arrayref of Bio::EnsEMBL::Feature (typically)
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

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
  return $self->merge_features(

    # use filter_by_min_max to filter out features we know are not overlapped
    # do it here so we're not affecting the cache
    $self->filter_features_by_min_max(\@features, @{$buffer->min_max})
  );
}


=head2 get_all_regions_by_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Example    : $regions = $as->get_all_regions_by_InputBuffer($ib)
  Description: Fetches all non-overlapping regions that overlap the variants currently
               in the input buffer. Regions are a simple arrayref [$chr, $start],
               corresponding to bins of genomic coordinates, size cache_region_size
  Returntype : arrayref of region arrayrefs
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

sub get_all_regions_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;
  my $cache_region_size = shift || $self->{cache_region_size};

  my @regions = ();
  my %seen = ();
  my $up_down_size = defined($self->{up_down_size}) ? $self->{up_down_size} : $self->up_down_size();

  my ($min, $max) = (1e10, 0);

  foreach my $vf(@{$buffer->buffer}) {
    my $chr = $vf->{chr} || $vf->slice->seq_region_name;
    throw("ERROR: Cannot get chromosome from VariationFeature") unless $chr;

    # find actual chromosome name used by AnnotationSource
    my $source_chr = $self->get_source_chr_name($chr);

    my ($vf_s, $vf_e) = ($vf->{start}, $vf->{end});

    # allow for indels
    ($vf_s, $vf_e) = ($vf_e, $vf_s) if $vf_s > $vf_e;

    # log min/max
    $min = $vf_s if $vf_s < $min;
    $max = $vf_e if $vf_e > $max;

    # convert to region-size
    my ($r_s, $r_e) = map {int(($_ - 1) / $cache_region_size)} ($vf_s - $up_down_size, $vf_e + $up_down_size);

    # add all regions between r_s and r_e inclusive
    for(my $s = $r_s; $s <= $r_e; $s++) {
      my $key = join(':', ($source_chr, $s));
      next if $seen{$key};

      push @regions, [$source_chr, $s];
      $seen{$key} = 1;
    }
  }

  $buffer->min_max([$min, $max]);

  return \@regions;
}


=head2 get_features_by_regions_cached

  Arg 1      : arrayref $regions
  Example    : $features = $as->get_features_by_regions_cached($regions)
  Description: Fetches features for given regions from in-memory cache
  Returntype : arrayref of Bio::EnsEMBL::Feature (typically)
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

sub get_features_by_regions_cached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;

  return [map {@{$cache->{$_->[0]}->{$_->[1]} || []}} @$regions];
}


=head2 filter_features_by_min_max

  Arg 1      : arrayref $features
  Arg 2      : int $min_min
  Arg 3      : int $max_coord
  Example    : $features = $as->filter_features_by_min_max($features, $min, $max)
  Description: Filters a list of features to those overlapping a given coordinate range.
               Used after features are fetched from regions/bins to limit those to be
               analysed to just those that overlap the min/max retrieved from the
               input buffer.
  Returntype : arrayref of Bio::EnsEMBL::Feature (typically)
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

sub filter_features_by_min_max {
  my $self = shift;
  my $features = shift;
  my $min = shift;
  my $max = shift;

  my $up_down_size = defined($self->{up_down_size}) ? $self->{up_down_size} : $self->up_down_size();
  $min -= $up_down_size;
  $max += $up_down_size;
  
  return [
    grep {overlap($_->{start}, $_->{end}, $min, $max)}
    @$features
  ];
}


=head2 filter_set

  Arg 1      : (optional) Bio::EnsEMBL::VEP::FilterSet $fs
  Example    : $fs = $as->filter_set()
  Description: Getter/setter for the FilterSet used to filter features
               retrieved by this annotation source.
  Returntype : Bio::EnsEMBL::VEP::FilterSet
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub filter_set {
  my $self = shift;
  $self->{filter_set} = shift if @_;
  return $self->{filter_set};
}


=head2 cache

  Example    : $cache = $as->cache()
  Description: Gets (and initialises if required) the cache hashref used
               to cache features between input buffer fills
  Returntype : hashref
  Exceptions : none
  Caller     : get_features_by_regions_cached()
  Status     : Stable

=cut

sub cache {
  return $_[0]->{_cache} ||= {};
}


=head2 clean_cache

  Example    : $as->clean_cache()
  Description: Removes out-of-range features from in-memory cache
  Returntype : none
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

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


=head2 info

  Example    : $info = $as->info()
  Description: Gets the info hashref for this annotation source, typically
               contains version information etc.
  Returntype : hashref
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::BaseRunner
  Status     : Stable

=cut

sub info {
  return $_[0]->{info} || {};
}

1;