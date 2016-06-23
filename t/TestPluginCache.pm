package TestPluginCache;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->{has_cache} = 1;

  return $self;
}

sub run {
  my ($self, $tva) = @_;
  $self->{cache}->{$tva->variation_feature->start} = 1;
  return { foo => 'bar' };
}

1;