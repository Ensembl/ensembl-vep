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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat - local disk regulatory feature annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat->new({
  config => $config,
  dir    => $dir
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Cache-based annotation source for regulatory data.

Data are stored as serialized objects on disk. Structure:

$cache = {
  chr => {
    RegulatoryFeature => [
      rf1,
      rf2
    ],
    MotifFeature => [
      mf1,
      mf2
    ]
  }
}

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::MotifFeature;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseSerialized
  Bio::EnsEMBL::VEP::AnnotationType::RegFeat
);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config $config,
                 dir    => string $dir,
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat
  Exceptions : none
  Caller     : CacheDir
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
  Description: Gets cell types available in this cache. Pre-loaded
               from cache info file.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : check_cell_types()
  Status     : Stable

=cut

sub get_available_cell_types {
  return $_[0]->{available_cell_types} || [];
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
    "%s/%s/%s_reg\.%s",
    $self->dir,
    $chr,
    $region,
    $self->file_suffix
  );
}


=head2 deserialized_obj_to_features

  Arg 1      : hashref $obj
  Example    : $features = $as->deserialized_obj_to_features($obj);
  Description: Takes the deserialized object read from the cache file,
               restores the objects by reattaching necessary adaptors.
  Returntype : arrayref of Bio::EnsEMBL::Feature
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub deserialized_obj_to_features {
  my $self = shift;
  my $obj = shift;

  my @features;
  my $sa = $self->get_adaptor('core', 'Slice');

  # the regfeat cache is structured like this:
  # $cache = {
  #   chr => {
  #     RegulatoryFeature => [
  #       rf1,
  #       rf2
  #     ],
  #     MotifFeature => [
  #       mf1,
  #       mf2
  #     ]
  #   }
  # }

  foreach my $chr(keys %$obj) {
    my $first = 1;

    # use reverse sort so RegulatoryFeature comes before MotifFeature
    foreach my $type(reverse sort keys %{$obj->{$chr}}) {
      foreach my $f(@{$obj->{$chr}->{$type}}) {

        # reattach slice adaptor
        # only need to do this for the first one as 
        if($first) {
          $f->{slice}->{adaptor} = $sa;
          $first = 0;
        }

        $f->{_vep_feature_type} ||= $type;

        push @features, $f;
      }
    }
  }

  return \@features;
}

1;