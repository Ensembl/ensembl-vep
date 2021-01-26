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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache - local disk annotation source base class.

=head1 SYNOPSIS

Should not be invoked directly

=head1 DESCRIPTION

Base class for annotation sources that read serialized data from a cache.

=head1 METHODS

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


=head2 deserialize_from_file

  Arg 1      : string $filename
  Example    : $obj = $as->deserialize_from_file($file);
  Description: Deserialize data from disk to in-memory object.
               Wrapper for deserialize_from_file_storable() or
               deserialize_from_file_sereal() depending on
               serializer used for cache.
  Returntype : scalar (typically hashref)
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub deserialize_from_file {
  my $self = shift;
  my $method = 'deserialize_from_file_'.$self->serializer_type;
  return $self->$method(@_);
}


=head2 deserialize_from_file_storable

  Arg 1      : string $filename
  Example    : $obj = $as->deserialize_from_file_storable($file);
  Description: Deserialize Storable data from disk to in-memory object.
               On-disk files are compressed with gzip, so this method
               will use the following methods of decompression in order
               of preference:
                - PerlIO::gzip
                - gzip binary
                - Compress::Zlib
  Returntype : scalar (typically hashref)
  Exceptions : none
  Caller     : deserialize_from_file()
  Status     : Stable

=cut

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


=head2 deserialize_from_file_sereal

  Arg 1      : string $filename
  Example    : $obj = $as->deserialize_from_file_sereal($file);
  Description: Deserialize Sereal data from disk to in-memory object.
  Returntype : scalar (typically hashref)
  Exceptions : none
  Caller     : deserialize_from_file()
  Status     : Stable

=cut

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


=head2 serializer_type

  Arg 1      : string $type
  Example    : $type = $as->serializer_type();
  Description: Getter/setter for serializer type - one of "storable" or "sereal"
  Returntype : string
  Exceptions : none
  Caller     : deserialize_from_file(), file_suffix()
  Status     : Stable

=cut

sub serializer_type {
  my $self = shift;
  
  $self->{serializer_type} = shift if @_;
  $self->{serializer_type} = 'storable' unless $self->{serializer_type};

  return $self->{serializer_type};
}


=head2 file_suffix

  Example    : $suffix = $as->file_suffix();
  Description: Get cache file suffix; depends on serializer_type
  Returntype : string
  Exceptions : none
  Caller     : get_dump_file_name()
  Status     : Stable

=cut

sub file_suffix {
  my $self = shift;
  return $self->{file_suffix} ||= $self->serializer_type eq 'sereal' ? 'sereal' : 'gz';
}


=head2 get_features_by_regions_uncached

  Arg 1      : arrayref $regions
  Example    : $features = $as->get_features_by_regions_uncached($regions)
  Description: Gets all features overlapping the given set of regions. See
               Bio::EnsEMBL::VEP::AnnotationSource::get_all_regions_by_InputBuffer()
               for information about regions.
  Returntype : arrayref
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

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
