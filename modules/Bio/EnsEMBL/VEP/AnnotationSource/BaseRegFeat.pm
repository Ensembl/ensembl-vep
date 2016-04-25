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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::BaseRegFeat
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::BaseRegFeat - base class for regfeat annotation sources

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::BaseRegFeat;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);

our @REG_FEAT_TYPES = qw(
  RegulatoryFeature
  MotifFeature
);

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  foreach my $rf(@{$self->get_all_features_by_InputBuffer($buffer, $self->{cache_region_size})}) {
    my $type = $rf->{_vep_feature_type} ||= (split('::', ref($rf)))[-1];
    my $constructor = 'Bio::EnsEMBL::Variation::'.$type.'Variation';
    my $add_method  = 'add_'.$type.'Variation';

    foreach my $vf(@{$buffer->get_overlapping_vfs($rf->{start}, $rf->{end})}) {
      $vf->{slice} ||= $rf->{slice};

      if(ref($vf) eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature') {
        my $svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
          -feature                      => $rf,
          -structural_variation_feature => $vf,
          -no_transfer                  => 1
        );

        push @{$vf->{regulation_structural_variations}->{$type}}, $svo;
      }
      
      else {
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
}

sub up_down_size {
  return 0;
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