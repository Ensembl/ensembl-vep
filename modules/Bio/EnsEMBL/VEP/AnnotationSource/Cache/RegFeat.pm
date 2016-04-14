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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat - local disk regulatory feature annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::MotifFeature;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache);

our @REG_FEAT_TYPES = qw(
  RegulatoryFeature
  MotifFeature
);

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  foreach my $rf(@{$self->get_all_features_by_InputBuffer($buffer, $self->{cache_region_size})}) {
    my $fs = $rf->{start};
    my $fe = $rf->{end};

    my $type = $rf->{_vep_feature_type} ||= (split('::', ref($rf)))[-1];
    my $constructor = 'Bio::EnsEMBL::Variation::'.$type.'Variation';
    my $add_method  = 'add_'.$type.'Variation';

    foreach my $vf(grep { overlap($fs, $fe, $_->{start}, $_->{end}) } @{$buffer->buffer}) {
      $vf->{slice} ||= $rf->{slice};

      $vf->$add_method(
        $constructor->new(
          -variation_feature  => $vf,
          -feature            => $rf,
          -no_ref_check       => 1,
          -no_transfer        => 1
        )
      );
    }
  }
}

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

    foreach my $type(keys %{$obj->{$chr}}) {
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

sub merge_features {
  my $self = shift;
  my $features = shift;

  my $seen = {};
  my @return;

  foreach my $f(@$features) {
    my $dbID = $f->dbID;
    my $type = $f->{_vep_feature_type};

    # the regfeat cache contains two feature types
    # each feature type has their own dbID namespace
    next if $seen->{$type}->{$dbID};

    push @return, $f;

    $seen->{$type}->{$dbID} = 1;
  }

  return \@return;
}

1;