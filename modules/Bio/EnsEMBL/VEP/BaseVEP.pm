=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(get_time);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::VEP::Stats;
use FileHandle;

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

# get stats object
sub stats {
  my $self = shift;

  my $config = $self->config;

  if(!exists($config->{_stats})) {
    $config->{_stats} = Bio::EnsEMBL::VEP::Stats->new({config => $self->config});
  }

  return $config->{_stats};
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

  my $config = $self->config;

  if(@_) {
    $config->{_registry} = shift;
  }
  elsif(!exists($config->{_registry})) {

    my $reg = 'Bio::EnsEMBL::Registry';
    
    if($self->param('database') || ($self->param('cache') && !$self->param('offline'))) {

      # load DB options from registry file if given
      if(my $registry_file = $self->param('registry')) {
        $self->status_msg("Loading DB self from registry file ", $registry_file) if $self->param('verbose');
        
        $reg->load_all(
          $registry_file,
          $self->param('verbose'),
          undef,
          $self->param('no_slice_cache')
        );
      }

      # otherwise manually connect to DB server
      else {
        my $species;

        unless($self->param('is_multispecies')) {
          if($self->param('species') =~ /^[a-z]+\_[a-z]+/i) {
            $species = $self->param('species');
          }
        }

        $reg->load_registry_from_db(
          -host       => $self->param('host'),
          -user       => $self->param('user'),
          -pass       => $self->param('password'),
          -port       => $self->param('port'),
          -db_version => $self->param('db_version'),
          -species    => $species,
          -verbose    => $self->param('verbose'),
          -no_cache   => $self->param('no_slice_cache'),
        );
      }

      eval { $reg->set_reconnect_when_lost() };
    }

    $config->{_registry} = $reg;
  }

  # copy to raw_config for plugin backward compatibility
  $config->{_params}->{reg} ||= $config->{_registry};

  return $config->{_registry};
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

    if($self->param('database') || ($self->param('cache') && !$self->param('offline'))) {
      $ad = $self->registry->get_adaptor($self->species, $group, $type)
    }

    $ad ||= $self->_get_fake_adaptor($group, $type);

    $cache->{$group}->{$type} = $ad;
  }

  return $cache->{$group}->{$type};
}

sub _get_fake_adaptor {
  my ($self, $group, $type) = @_;

  my $ad;

  my $module_name = sprintf(
    "Bio::EnsEMBL::%sDBSQL::%sAdaptor",
    (lc($group) eq 'core' ? '' : ucfirst($group).'::'),
    $type
  );

  eval { $ad = $module_name->new_fake($self->species()) };

  return $ad;
}

sub get_slice {
  my $self = shift;
  my $chr = shift;

  my $cache = $self->{_slice_cache} ||= {};

  if(!exists($cache->{$chr})) {
    my $slice;

    if(my $sa = $self->get_adaptor($self->param('core_type'), 'Slice')) {
      $slice = $sa->fetch_by_region(undef, $chr);
    }

    elsif(my $fasta_db = $self->fasta_db) {
      my $fa_length = $fasta_db->length($chr);
      my $length = $fa_length && $fa_length > 0 ? $fa_length : 1;

      $slice = Bio::EnsEMBL::Slice->new(
        -COORD_SYSTEM      => $self->{_coord_system} ||= Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome', -RANK => 1),
        -START             => 1,
        -END               => $length,
        -SEQ_REGION_NAME   => $chr,
        -SEQ_REGION_LENGTH => $length
      );

      $slice->{is_fake} = 1;
    }

    $cache->{$chr} = $slice;
  }

  return $cache->{$chr};
}

sub get_database_assembly {
  my $self = shift;

  my $config = $self->config;

  if(!exists($config->{_database_assembly})) {
    my $assembly;

    if(my $csa = $self->get_adaptor('core', 'CoordSystem')) {
      my ($highest_cs) = @{$csa->fetch_all()};
      $assembly = $highest_cs->version();
    }

    $config->{_database_assembly} = $assembly;
  }

  return $config->{_database_assembly};
}

sub fasta_db {
  my $self = shift;

  if(!exists($self->config->{_fasta_db})) {
    my $fasta_db;

    if(my $fasta_file = $self->param('fasta')) {

      $fasta_db = setup_fasta(
        -FASTA => $fasta_file,
        -ASSEMBLY => $self->param('assembly'),
        -OFFLINE => $self->param('offline'),
        -SYNONYMS => $self->chromosome_synonyms,
      );

      $self->stats->log_fasta_chromosomes($fasta_db);
    }

    $self->config->{_fasta_db} = $fasta_db;
  }

  return $self->config->{_fasta_db};
}

sub chromosome_synonyms {
  my $self = shift;
  my $file = shift;

  if($file) {
    open IN, $file or throw("ERROR: Could not read synonyms file $file: $!");

    my $synonyms = $self->config->{_chromosome_synonyms} ||= {};

    while(<IN>) {
      chomp;
      my @split = split(/\s+/, $_);

      my $ref = shift @split;

      foreach my $syn(@split) {
        $synonyms->{$ref}->{$syn} = 1;
        $synonyms->{$syn}->{$ref} = 1;
      }
    }

    close IN;
  }

  return $self->config->{_chromosome_synonyms} ||= {};
}

# adds shortcuts to named params to this object
sub add_shortcuts {
  my $self = shift;

  return unless $self->config;

  foreach my $param(map {ref($_) eq 'ARRAY' ? @$_ : $_} @_) {
    throw("ERROR: add_shortcuts would overwrite value for \"$param\"\n") if exists($self->{$param});
    $self->{$param} = $self->param($param); 
  }
}

# prints a status message to STDOUT
sub status_msg {
  my $self = shift;
  return if $self->param('quiet');
  
  my $msg = (@_ ? (join "", @_) : "No message");
  print get_time()." - ".$msg.($msg =~ /\n$/ ? "" : "\n");
}

# prints warning messages to STDERR or a log file
sub warning_msg {
  my $self = shift;
  my $text = shift || '';

  $text = 'WARNING: '.$text unless $text =~ /^warn/i;
  $text = $text."\n" unless $text =~ /\n$/;

  # $self->config->{warning_count}++ if $self->config;

  my $fh = $self->warning_fh;

  print $fh $text;

  unless($self->param('quiet')) {
    warn($self->param('no_progress') ? $text : "\n$text");
  }
}

# gets filehandle for use by warning_msg
sub warning_fh {
  my $self = shift;

  unless(exists($self->config->{warning_fh})) {
    my $file = $self->param('warning_file');
    my $fh;

    if($file && $file =~ /^stderr$/i) {
      $fh = *STDERR;
    }

    else {
      $file ||= ($self->param('output_file') || 'vep').'_warnings.txt';
      $self->param('warning_file', $file);

      $fh = FileHandle->new();
      $fh->open(">".$file) or throw("ERROR: Could not write to warnings file $file\n");
    }

    $self->config->{warning_fh} = $fh;
  }

  return $self->config->{warning_fh};
}

sub module_prefix {
  my $self = shift;
  my $prefix = 'Bio::EnsEMBL::VEP';
  $prefix .= '::Haplo' if $self->param('haplo');
  return $prefix;
}

1;
