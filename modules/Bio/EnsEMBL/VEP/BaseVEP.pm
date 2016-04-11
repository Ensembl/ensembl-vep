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

# EnsEMBL module for Bio::EnsEMBL::VEP::BaseVEP
#
#

=head1 NAME

Bio::EnsEMBL::VEP::BaseVEP - Base class used by the Bio::EnsEMBL::VEP::* classes

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::BaseVEP;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# new method, may or may not be reused by child classes
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # initialise self
  my $self = bless {}, $class;
  
  # most new methods will pass in the config in the hashref given as the first arg to new()
  if(my $hashref = $_[0]) {

    # does the config key exist?
    if(exists($hashref->{config})) {

      # copy and delete it so child class new methods don't get confused
      # if they want to deal with this differently, then they shouldn't run $class->SUPER::new
      my $config = delete($hashref->{config});

      # do an assert ref to check the class
      assert_ref($config, 'Bio::EnsEMBL::VEP::Config');

      $self->{_config} = $config;
    }
  }

  return $self;
}

# returns Bio::EnsEMBL::Variation::VEP::Config object
sub config {
  return $_[0]->{_config};
}

# gets/sets the value of a config parameter given a key name
sub param {
  my $self = shift;
  return $self->config->param(@_);
}

sub species {
  my $self = shift;

  if(@_) {
    $self->{_species} = shift;
    $self->param('species', $self->{_species});
  }
  else {
    $self->{_species} ||= $self->param('species');
  }

  return $self->{_species};
}

# gets Bio::EnsEMBL::Registry, connects to DB if available
sub registry {
  my $self = shift;

  if(!exists($self->{_registry})) {

    my $reg = 'Bio::EnsEMBL::Registry';
    
    unless($self->param('offline')) {
      # suppress warnings that the FeatureAdpators spit if using no_slice_cache
      Bio::EnsEMBL::Utils::Exception::verbose(1999) if $self->param('no_slice_cache');

      # load DB options from registry file if given
      if(my $registry_file = $self->param('registry')) {
        $self->status_msg("Loading DB self from registry file ", $registry_file);
        
        $reg->load_all(
          $registry_file,
          $self->param('verbose'),
          undef,
          $self->param('no_slice_cache')
        );
      }

      # otherwise manually connect to DB server
      else {
        $reg->load_registry_from_db(
          -host       => $self->param('host'),
          -user       => $self->param('user'),
          -pass       => $self->param('password'),
          -port       => $self->param('port'),
          -db_version => $self->param('db_version'),
          -species    => $self->param('species') =~ /^[a-z]+\_[a-z]+/i ? $self->param('species') : undef,
          -verbose    => $self->param('verbose'),
          -no_cache   => $self->param('no_slice_cache'),
        );
      }

      eval { $reg->set_reconnect_when_lost() };
    }

    $self->{_registry} = $reg;
  }

  return $self->{_registry};
}

sub get_adaptor {
  my $self = shift;
  my $group = shift;
  my $type = shift;

  throw("No adaptor group specified") unless $group;
  throw("No adaptor type specified") unless $type;

  my $cache = $self->{_adaptors} ||= {};

  if(!exists($cache->{$group}) || !exists($cache->{$group}->{$type})) {
    my $ad;

    if($self->param('offline')) {
      my $module_name = sprintf(
        "Bio::EnsEMBL::%sDBSQL::%sAdaptor",
        (lc($group) eq 'core' ? '' : ucfirst($group).'::'),
        $type
      );

      eval { $ad = $module_name->new_fake($self->species()) };
    }

    else {
      $ad = $self->registry->get_adaptor($self->species, $group, $type);
    }

    $cache->{$group}->{$type} = $ad;
  }

  return $cache->{$group}->{$type};
}

# prints a status message to STDOUT
sub status_msg {
  my $self = shift;
  my $config = $self->config();
  
  return if defined($config->{quiet});
  
  my $msg = (@_ ? (join "", @_) : "No message");
  print $self->get_time()." - ".$msg.($msg =~ /\n$/ ? "" : "\n");
}

# gets time
sub get_time() {
  my @time = localtime(time());

  # increment the month (Jan = 0)
  $time[4]++;

  # add leading zeroes as required
  for my $i(0..4) {
    $time[$i] = "0".$time[$i] if $time[$i] < 10;
  }

  # put the components together in a string
  my $time =
    ($time[5] + 1900)."-".
    $time[4]."-".
    $time[3]." ".
    $time[2].":".
    $time[1].":".
    $time[0];

  return $time;
}

1;
