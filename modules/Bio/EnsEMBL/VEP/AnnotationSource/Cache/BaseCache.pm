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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache - local disk annotation source base class.

DO NOT USE DIRECTLY

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache;

use Storable qw(nstore_fd fd_retrieve freeze thaw);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::Cache);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(compress cache_region_size)]);

  return $self;
}

sub deserialize_from_file {
  my $self = shift;
  
  my $type = $self->serializer_type;
  my $method = 'deserialize_from_file_'.$type;

  return $self->$method(@_);
}

sub deserialize_from_file_storable {
  my $self = shift;
  my $file = shift;

  open my $fh, $self->{compress}." ".$file." |" or die "ERROR: $!";
  my $obj = fd_retrieve($fh);
  close $fh;

  return $obj;
}

sub deserialize_from_file_sereal {
  my $self = shift;
  my $file = shift;

  $self->{decoder} ||= Sereal::Decoder->new();
  open IN, $file;
  my $obj = $self->{decoder}->decode(join('', <IN>));
  close IN;

  return $obj;
}

sub serializer_type {
  my $self = shift;
  
  $self->{serializer_type} = shift if @_;
  $self->{serializer_type} = 'storable' unless $self->{serializer_type};

  return $self->{serializer_type};
}

sub file_suffix {
  my $self = shift;
  return $self->{file_suffix} ||= $self->serializer_type eq 'sereal' ? 'sereal' : 'gz';
}

sub get_all_features_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my $regions = $buffer->get_cache_regions($self->{cache_region_size});

  my @not_cached;
  my @features;

  # attempt to fetch from memory first
  foreach my $region(@$regions) {
    my $cache_features = $self->get_features_by_regions_from_memory([$region]);

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
  push @features, @{$self->get_features_by_regions_from_disk(\@not_cached)};

  $self->clean_cache($regions);

  # use merge_features to remove duplicates
  return $self->merge_features(\@features);
}

sub get_features_by_regions_from_memory {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;

  return [map {@{$cache->{$_->[0]}->{$_->[1]} || []}} @$regions];
}

sub get_features_by_regions_from_disk {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;
  my $cache_region_size = $self->{cache_region_size};
  my @features;

  foreach my $region(@{$regions}) {
    my ($c, $s) = @$region;

    my $file = $self->get_dump_file_name(
      $c,
      ($s * $cache_region_size) + 1,
      ($s + 1) * $cache_region_size
    );

    push @features, @{$self->deserialized_obj_to_features(
      $self->deserialize_from_file($file)
    )} if -e $file;

    $cache->{$c}->{$s} = \@features;
  }

  return \@features;
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
