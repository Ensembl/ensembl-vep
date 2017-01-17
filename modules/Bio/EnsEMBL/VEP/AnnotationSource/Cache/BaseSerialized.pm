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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache - local disk annotation source base class.

DO NOT USE DIRECTLY

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseSerialized;

use Storable qw(nstore_fd fd_retrieve);
use Compress::Zlib;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::Cache);

our ($CAN_USE_PERLIO_GZIP, $CAN_USE_GZIP, $CAN_USE_SEREAL);

BEGIN {

  # check Sereal
  if (eval q{ require Sereal; 1 }) {
    $CAN_USE_SEREAL = 1;
  }

  # check PerlIO::gzip
  if (eval q{ require PerlIO::gzip; 1 }) {
    $CAN_USE_PERLIO_GZIP = 1;
  }

  # check gzip
  if (`which gzip` =~ /\/gzip/) {
    $CAN_USE_GZIP = 1;
  }
}

sub deserialize_from_file {
  my $self = shift;
  my $method = 'deserialize_from_file_'.$self->serializer_type;
  return $self->$method(@_);
}

sub deserialize_from_file_storable {
  my $self = shift;
  my $file = shift;

  # we have three options to decompress, try in order of speed:
  my $obj;

  # 1) PerlIO::gzip
  if($CAN_USE_PERLIO_GZIP) {
    open my $fh, "<:gzip", $file or throw("ERROR: $!");
    $obj = fd_retrieve($fh);
  }

  # 2) gzip binary
  elsif($CAN_USE_GZIP) {
    open my $fh, "gzip -dc $file |" or throw("ERROR: $!");
    $obj = fd_retrieve($fh);
  }

  # 3) Compress::Zlib
  else {
    my $gz = gzopen($file, 'rb') or throw("ERROR: $!");

    my ($buffer, $serialized);
    while($gz->gzread($buffer)) {
      $serialized .= $buffer;
    }

    # now use fd_retrieve on a made-up filehandle to deserialize
    open IN, '<', \$serialized;
    $obj = fd_retrieve(\*IN);
    close IN;
  }

  return $obj;
}

sub deserialize_from_file_sereal {
  my $self = shift;
  my $file = shift;

  throw("ERROR: Unable to use Sereal; module not installed\n") unless $CAN_USE_SEREAL;

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

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;
  my $cache_region_size = $self->{cache_region_size};

  my @return;

  foreach my $region(@{$regions}) {
    my ($c, $s) = @$region;

    my $file = $self->get_dump_file_name(
      $c,
      ($s * $cache_region_size) + 1,
      ($s + 1) * $cache_region_size
    );

    my @features = @{$self->deserialized_obj_to_features(
      $self->deserialize_from_file($file)
    )} if -e $file;

    $cache->{$c}->{$s} = \@features;

    push @return, @features;
  }

  return \@return;
}

1;
