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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationType::RegFeat
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationType::RegFeat - base class for regfeat annotation sources

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Helper class for all RegFeat-based AnnotationSource classes. Contains the bulk of the
code that carries out the annotation.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationType::RegFeat;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

our @REG_FEAT_TYPES = qw(
  RegulatoryFeature
  MotifFeature
);


=head2 check_cell_types

  Example    : $ok = $as->check_cell_types();
  Description: Check if this annotation source contains the cell types
               as specified by the user with --cell_type.
  Returntype : bool
  Exceptions : throws if cell type unavailable
  Caller     : new()
  Status     : Stable

=cut

sub check_cell_types {
  my $self = shift;

  my $given_cell_types = $self->{cell_type};

  return unless $given_cell_types && scalar @$given_cell_types;

  my %available = map {$_ => 1} @{$self->get_available_cell_types};

  foreach my $ct(@$given_cell_types) {
    throw(
      sprintf(
        "ERROR: Cell type \"%s\" unavailable; available cell types are:\n%s\n",
        $ct,
        join(",", @{$self->get_available_cell_types})
      )
    ) unless $available{$ct};
  }

  return 1;
}


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $as->annotate_InputBuffer($ib);
  Description: Gets overlapping regulatory and motif features for the variants in
               the input buffer, and creates VariationFeatureOverlap objects for
               each overlapping pair of variant/feature.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  foreach my $rf(@{$self->get_all_features_by_InputBuffer($buffer, $self->{cache_region_size})}) {
    my $type = $rf->{_vep_feature_type} ||= (split('::', ref($rf)))[-1];
    my $constructor = 'Bio::EnsEMBL::Variation::'.$type.'Variation';
    my $add_method  = 'add_'.$type.'Variation';

    my $slice = $rf->{slice};

    foreach my $vf(@{$buffer->get_overlapping_vfs($rf->{start}, $rf->{end})}) {
      $vf->{slice} ||= $slice;

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


=head2 up_down_size

  Example    : $size = $as->up_down_size();
  Description: Gets range in bp that should be added to boundaries
               when fetching features. 0 for regulatory features.
  Returntype : int
  Exceptions : none
  Caller     : get_all_regions_by_InputBuffer()
  Status     : Stable

=cut

sub up_down_size {
  return 0;
}


=head2 merge_features

  Arg 1      : arrayref of Bio::EnsEMBL::Feature $features
  Example    : $merged = $as->merge_features($features);
  Description: "Uniquifies" a list of features using dbID
  Returntype : arrayref of Bio::EnsEMBL::Feature
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

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