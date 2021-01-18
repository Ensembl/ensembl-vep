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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat - database RegFeat annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat->new({
  config => $config,
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Database-based annotation source for regulatory data.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat;

use Scalar::Util qw(weaken isweak);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Database
  Bio::EnsEMBL::VEP::AnnotationType::RegFeat
);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config $config,
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat
  Exceptions : none
  Caller     : AnnotationSourceAdaptor
  Status     : Stable

=cut

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


=head2 get_available_cell_types

  Example    : $types = $as->get_available_cell_types();
  Description: Gets cell types available in this database.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : check_cell_types()
  Status     : Stable

=cut

sub get_available_cell_types {
  my $self = shift;

  if(!exists($self->{available_cell_types})) {
    my $regulatory_build_adaptor = $self->get_adaptor('funcgen', 'RegulatoryBuild');
    my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
    $self->{available_cell_types} = [
      sort
      map {s/ /\_/g; $_}
      map {$_->short_name}
      @{$regulatory_build->get_all_Epigenomes}
    ];
  }

  return $self->{available_cell_types};
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

    $sub_slice->{coord_system}->{adaptor} ||= $self->get_adaptor('core', 'coordsystem');
    $sub_slice->{adaptor} ||= $self->get_adaptor('core', 'slice');
    my $type = 'RegulatoryFeature'; 
    my $features = $self->get_adaptor('funcgen', $type)->fetch_all_by_Slice($sub_slice);

    next unless defined($features);

    foreach my $rf(@$features) { 
      # weaken a circular ref
      foreach my $a(@{$rf->{_regulatory_activity} || []}) {
        weaken($a->{_regulatory_feature}) unless isweak($a->{_regulatory_feature});
      }

      # cell types
      if($self->{cell_type} && scalar(@{$self->{cell_type}})) {
        my %cl =
          map {$_->[0] => $_->[1]}
          map {$_->[0] =~ s/ /\_/g; $_}
          map {[$_->get_Epigenome->short_name, $_->activity]}
          @{$rf->regulatory_activity};
        $rf->{cell_types} = \%cl;
      }
    }
    push @region_features,
      map { $_->{_vep_feature_type} ||= $type; $_ }
      map { $_->transfer($sr_slice) }
      @{$features};

    $type = 'MotifFeature';
    my @motif_features = ();
    foreach my $rf (@$features) {
      my $rf_seq_region = $rf->seq_region_name;
      my @experimentally_verified_MotifFeatures = grep {$_->seq_region_name eq "$rf_seq_region"} @{$rf->get_all_experimentally_verified_MotifFeatures()};  
      push @motif_features, @experimentally_verified_MotifFeatures;
    } 

    foreach my $mf(@motif_features) {
      $mf->get_BindingMatrix->summary_as_hash();
      if($self->{cell_type} && scalar(@{$self->{cell_type}})) {
        my %cl =  
          map {$_->[0] => $_->[1]}
          map {$_->[0] =~ s/ /\_/g; $_}
          map { [$_->get_PeakCalling->get_Epigenome->name, 1] } @{$mf->get_all_overlapping_Peaks};
        $mf->{cell_types} = \%cl;
      }
    }
    push @region_features,
      map { $_->{_vep_feature_type} ||= $type; $_ }
      map { $_->transfer($sr_slice) }
      @motif_features;

    $cache->{$c}->{$region_start} = \@region_features;
    push @return, @region_features;
  }

  return \@return;
}


=head2 info

  Example    : $info = $as->info()
  Description: Gets the info hashref for this annotation source. Contains
               the regulatory_build version.
  Returntype : hashref
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::BaseRunner
  Status     : Stable

=cut

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
