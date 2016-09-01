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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation - database structural variation annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation;

use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::Database);

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  foreach my $svf(@{$self->get_all_features_by_InputBuffer($buffer)}) {
    foreach my $vf(@{$buffer->get_overlapping_vfs($svf->{start}, $svf->{end})}) {
      $vf->{overlapping_svs}->{$svf->variation_name} = 1;
    }
  }
}

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;
  my @return;

  my $cache_region_size = $self->{cache_region_size};

  my $svfa = $self->get_adaptor('variation', 'StructuralVariationFeature');
  my $sa = $self->get_adaptor('core', 'Slice');

  foreach my $region(@$regions) {

    my ($c, $region_start) = @{$region};

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

    push @return,
      map {$_->transfer($sr_slice)}
      @{$svfa->fetch_all_by_Slice($sub_slice)};
  }

  return \@return;
}

sub merge_features {
  my $self = shift;
  my $features = shift;

  my %seen = ();
  my @return;

  foreach my $feature(@$features) {
    push @return, $feature unless $seen{$feature->dbID};
    $seen{$feature->dbID} = 1;
  }

  return \@return;
}

sub up_down_size {
  return 1;
}

1;