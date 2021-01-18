=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Regulation;

use strict;
use warnings;

use Scalar::Util qw(weaken);
use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat;
use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper);

sub run {
  my $self = shift;

  my $vep_params = $self->get_vep_params();

  my $config = Bio::EnsEMBL::VEP::Config->new($vep_params);

  my $region_size = $self->param('region_size');
  my $hive_dbc = $self->dbc;
  $hive_dbc->disconnect_if_idle() if defined $hive_dbc;

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat->new({
    config => $config,
    cache_region_size => $region_size,
  });

  # manually turn on cell type fetch, we have to do it here to avoid the check killing it
  $as->{cell_type} = [1];

  my $cache = Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat->new({
    config => $config,
    cache_region_size => $region_size,
    dir => $self->get_cache_dir($vep_params)
  });

  $self->dump_chrs($as, $cache);

  $self->dump_info($as, $self->get_cache_dir($vep_params));
  
  return;
}

sub dump_info {
  my ($self, $as, $dir) = @_;

  my $info_file = $dir.'/info.txt_regulation';
  return if -e $info_file;

  open OUT, ">$info_file";

  # indicate contains regulation data
  print OUT "regulatory\t1\n";

  # cell types
  my $cell_types = $as->get_available_cell_types;
  if($cell_types && @$cell_types) {
    print OUT "cell_types\t".join(",", @$cell_types)."\n";
  }

  my $info = $as->info;
  print OUT "source_$_\t".$info->{$_}."\n" for keys %$info;

  close OUT;
}

sub get_dumpable_object {
  my ($self, $as, $sr, $chr, $s) = @_;

  my %split;
  push @{$split{$_->{_vep_feature_type}}}, $_ for
    map {$self->clean_regfeat($_)}
    @{$as->get_features_by_regions_uncached([[$sr, $s]], 1)};

  my $obj = { $chr => \%split };

  # delete adaptors here before dumping
  # they get reattached to the slice within get_features_by_regions_uncached
  foreach my $type(keys %{$obj->{$chr}}) {
    if(my $rf = $obj->{$chr}->{$type}->[0]) {
      delete $rf->{slice}->{adaptor};
      delete $rf->{slice}->{coord_system}->{adaptor};
    }
  }

  return $obj;
}

sub clean_regfeat {
  my $self = shift;
  my $rf = shift;

  delete $rf->{$_} for qw(
    adaptor
    binary_string
    bound_start
    bound_end
    attribute_cache
    feature_set
    analysis
    set
    _regulatory_activity
    _regulatory_build
    overlapping_Peaks
  );

  if(defined($rf->{binding_matrix})) {
    $rf->{_variation_effect_feature_cache}->{seq} = $rf->seq;
    foreach my $key(qw(adaptor feature_type analysis dbID)) {
      delete $rf->{binding_matrix}->{$key};
    }
    if (defined($rf->{binding_matrix}->{associated_transcription_factor_complexes})) {
      foreach my $tfc (@{$rf->{binding_matrix}->{associated_transcription_factor_complexes}}) {
        foreach my $key(qw(adaptor dbID feature_type)) {
          delete $tfc->{$key};
        }
        foreach my $component (@{$tfc->{components}}) {
          foreach my $key(qw(adaptor dbID feature_type)) {
            delete $component->{$key};
          }
        }
      }
    }
  }

  $rf->{feature_type} = $rf->feature_so_term if $rf->{feature_type};

  return $rf;
}

1;
