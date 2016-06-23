package TestPluginRegFeat;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
  return ['RegulatoryFeature'];
}

1;