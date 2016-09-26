=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::VEP::FilterSet;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

sub new {  
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  $self->{$_} = $hashref->{$_} for keys %$hashref;

  $self->filter_set(Bio::EnsEMBL::VEP::FilterSet->new($hashref->{filter})) if $hashref->{filter};

  return $self;
}

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
    my @region_starts;

    if($vf_s > $vf_e) {
      @region_starts = ($vf_e - $up_down_size, $vf_s + $up_down_size);

      $min = $vf_e if $vf_e < $min;
      $max = $vf_s if $vf_s > $max;
    }
    else {
      @region_starts = ($vf_s - $up_down_size, $vf_e + $up_down_size);

      $min = $vf_s if $vf_s < $min;
      $max = $vf_e if $vf_e > $max;
    }

    foreach my $region_start(map {int(($_ - 1)/ $cache_region_size)} @region_starts) {
      my $key = join(':', ($source_chr, $region_start));
      next if $seen{$key};

      push @regions, [$source_chr, $region_start];
      $seen{$key} = 1;
    }
  }

  $buffer->min_max([$min, $max]);

  return \@regions;
}

sub get_source_chr_name {
  my $self = shift;
  my $chr = shift;

  my $chr_name_map = $self->{_chr_name_map} ||= {};

  if(!exists($chr_name_map->{$chr})) {
    my $mapped_name = $chr;

    my %valid = map {$_ => 1} @{$self->can('get_valid_chromosomes') ? $self->get_valid_chromosomes : []};

    unless($valid{$chr}) {

      # try synonyms first
      my $synonyms = $self->chromosome_synonyms;

      foreach my $syn(keys %{$synonyms->{$chr} || {}}) {
        if($valid{$syn}) {
          $mapped_name = $syn;
          last;
        }
      }

      # still haven't got it
      if($mapped_name eq $chr) {

        # try adding/removing "chr"
        if($chr =~ /^chr/i) {
          my $tmp = $chr;
          $tmp =~ s/^chr//i;

          $mapped_name = $tmp if $valid{$tmp};
        }
        elsif($valid{'chr'.$chr}) {
          $mapped_name = 'chr'.$chr;
        }
      }
    }

    $chr_name_map->{$chr} = $mapped_name;
  }

  return $chr_name_map->{$chr};
}

sub get_features_by_regions_cached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;

  return [map {@{$cache->{$_->[0]}->{$_->[1]} || []}} @$regions];
}

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

sub filter_set {
  my $self = shift;
  $self->{filter_set} = shift if @_;
  return $self->{filter_set};
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

sub info {
  return $_[0]->{info} || {};
}

1;