package TestPlugin;
use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub run {
  if($_[-1] eq 'skip') {
    return undef;
  }
  elsif($_[-1] eq 'not_hash') {
    return [ 'Hello' ];
  }
  else {
    return { test => 'Hello' };
  }
}

sub get_header_info {
  return { test => 'header' };
}

1;