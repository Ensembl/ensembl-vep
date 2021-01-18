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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Variation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Variation - class for cache variation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation->new({
  config => $config,
  dir    => $dir
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Cache-based annotation source for known variant data. Variants are
stored in 1MB chunks on disk, one line per variant.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation;

use Scalar::Util qw(weaken);
use Compress::Zlib;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation
);


=head2 get_features_by_regions_uncached

  Arg 1      : arrayref $regions
  Example    : $vars = $as->get_features_by_regions_uncached($regions)
  Description: Gets all known variants overlapping the given set of regions. See
               Bio::EnsEMBL::VEP::AnnotationSource::get_all_regions_by_InputBuffer()
               for information about regions.
  Returntype : arrayref of variant hashes
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

    my @features = @{$self->read_variations_from_file($file)} if -e $file;

    $cache->{$c}->{$s} = \@features;

    push @return, @features;
  }

  return \@return;
}


=head2 read_variations_from_file

  Arg 1      : string $file
  Example    : $vars = $as->get_features_by_regions_uncached($regions)
  Description: Reads a cache file of variants, creating a list of variant
               hash records
  Returntype : arrayref of variant hashes
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub read_variations_from_file {
  my $self = shift;
  my $file = shift;

  my @vars;

  # use Compress::Zlib interface
  my $gz = gzopen($file, 'rb');

  my $line;
  while($gz->gzreadline($line)) {
    chomp($line);
    my $var = $self->parse_variation($line);
    push @vars, $var if $self->filter_variation($var);
  }

  return \@vars;
}


=head2 get_dump_file_name

  Arg 1      : string $chr
  Arg 2      : string $region or int $start
  Arg 3      : (optional) int $end
  Example    : $file = $as->get_dump_file_name(1, "1-1000000");
  Description: Gets file name from the cache given a region.
  Returntype : string
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub get_dump_file_name {
  my $self = shift;
  my $chr  = shift;
  my $region = shift;

  throw("No chromosome given") unless $chr;
  throw("No region given") unless $region;

  # allow to pass region (start-end) or $start, $end
  $region .= '-'.shift if @_;

  return sprintf(
    "%s/%s/%s_var\.gz",
    $self->dir,
    $chr,
    $region,
  );
}

1;
