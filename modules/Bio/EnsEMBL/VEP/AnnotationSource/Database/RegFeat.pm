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

  $self->check_cell_types();

  return $self;
}

sub get_available_cell_types {
  my $self = shift;

  if(!exists($self->{available_cell_types})) {
    my $regulatory_build_adaptor = $self->get_adaptor('funcgen', 'RegulatoryBuild');
    my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
    $self->{available_cell_types} = [
      sort
      map {s/ /\_/g; $_}
      map {$_->display_label}
      @{$regulatory_build->get_all_Epigenomes}
    ];
  }

  return $self->{available_cell_types};
}

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;
  my $chr_is_seq_region = shift;

  my $cache = $self->cache;
  my @return;

  my $cache_region_size = $self->{cache_region_size};

  foreach my $region(@{$regions}) {
    my ($c, $region_start) = @$region;

    my $slice = $self->get_slice($c, $chr_is_seq_region);

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

    my $rfa = $self->get_adaptor('funcgen', 'RegulatoryFeature');

    foreach my $type(qw(RegulatoryFeature MotifFeature)) {
      my $features = $self->get_adaptor('funcgen', $type)->fetch_all_by_Slice($sub_slice);
      next unless defined($features);

      # cell types
      if($self->{cell_type} && scalar(@{$self->{cell_type}})) {
        foreach my $rf(@$features) {

          my %cl;

          # get cell type using regulatory_activity objects
          if($type eq 'RegulatoryFeature') {
            %cl =
              map {$_->[0] => $_->[1]}
              map {$_->[0] =~ s/ /\_/g; $_}
              map {[$_->epigenome->display_label, $_->activity]}
              grep {!$_->_is_multicell}
              @{$rf->regulatory_activity};
          }

          # get cell type by fetching regfeats that contain this MotifFeature
          elsif($type eq 'MotifFeature') {
            %cl =
              map {$_->[0] => $_->[1]}
              map {$_->[0] =~ s/ /\_/g; $_}
              map {[$_->epigenome->display_label, $_->activity]}
              grep {!$_->_is_multicell}
              map {@{$_->regulatory_activity}}
              @{$rfa->fetch_all_by_attribute_feature($rf)};
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

sub info {
  my $self = shift;

  if(!exists($self->{info})) {
    my $version;

    if(my $fg_mca = $self->get_adaptor('funcgen', 'metacontainer')) {
      foreach my $meta_key(qw(regbuild.version)) {
        my $meta_version = $fg_mca->list_value_by_key($meta_key);
        $version = $meta_version->[0] if defined($meta_version) && scalar @$meta_version;
      }

      # from 85 version is in regulatory_build table
      unless($version) {
        my $sth = $fg_mca->db->dbc->prepare(qq{
          SELECT version FROM regulatory_build 
        });
        $sth->execute;

        $sth->bind_columns(\$version);
        $sth->fetch;
        $sth->finish;
      }
    }

    $self->{info}->{regbuild} = $version;
  }

  return $self->{info};
}

1;