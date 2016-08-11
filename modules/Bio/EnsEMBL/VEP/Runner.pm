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

# EnsEMBL module for Bio::EnsEMBL::VEP::Runner
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Runner - runner class for VEP

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Runner;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Storable qw(freeze thaw);
use IO::Socket;
use IO::Select;
use FileHandle;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(get_time merge_hashes);
use Bio::EnsEMBL::VEP::Constants;
use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::Parser;
use Bio::EnsEMBL::VEP::InputBuffer;
use Bio::EnsEMBL::VEP::OutputFactory;
use Bio::EnsEMBL::VEP::AnnotationSourceAdaptor;

# don't assert refs
$Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0;

# don't use rearrange
$Bio::EnsEMBL::Utils::Argument::NO_REARRANGE = 1;

# avoid using transfer
$Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER = 1;

# has our own new method, does not use BaseVEP's
# since this is the class users will be instantiating
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # initialise self
  my $self = bless {}, $class;

  # get a config object
  $self->{_config} = Bio::EnsEMBL::VEP::Config->new(@_);

  return $self;
}

# dispatcher/runner for all initial setup from config
sub init {
  my $self = shift;

  return 1 if $self->{_initialized};

  # log start time
  $self->stats->start_time();

  # setup DB connection
  $self->setup_db_connection();

  # get chromosome synoyms
  $self->chromosome_synonyms($self->param('synonyms'));

  my $plugins = $self->get_all_Plugins();

  # get all annotation sources
  my $annotation_sources = $self->get_all_AnnotationSources();

  # setup FASTA file DB
  $self->fasta_db();

  my $buffer = $self->get_InputBuffer();

  $self->post_setup_checks();

  $self->stats->info($self->get_output_header_info);

  return $self->{_initialized} = 1;
}

# run
sub run {
  my $self = shift;

  $self->init();

  my $fh = $self->get_output_file_handle();

  print $fh "$_\n" for @{$self->get_OutputFactory->headers};

  while(my $line = $self->next_output_line) {
    print $fh "$line\n";
  }

  close $fh;

  $self->dump_stats;

  return 1;
}

sub next_output_line {
  my $self = shift;

  my $output_buffer = $self->{_output_buffer} ||= [];

  return shift @$output_buffer if @$output_buffer;

  $self->init();

  if($self->param('fork')) {
    push @$output_buffer, @{$self->_forked_buffer_to_output($self->get_InputBuffer)};
  }
  else {
    push @$output_buffer, @{$self->_buffer_to_output($self->get_InputBuffer)};
  }

  return @$output_buffer ? shift @$output_buffer : undef;
}

sub _buffer_to_output {
  my $self = shift;
  my $input_buffer = shift;

  my @output;
  my $vfs = $input_buffer->next();

  if($vfs && scalar @$vfs) {
    my $output_factory = $self->get_OutputFactory;

    foreach my $as(@{$self->get_all_AnnotationSources}) {
      $as->annotate_InputBuffer($input_buffer);
    }
      
    $input_buffer->finish_annotation;

    push @output, @{$output_factory->get_all_lines_by_InputBuffer($input_buffer)};
  }

  return \@output;
}

sub _forked_buffer_to_output {
  my $self = shift;
  my $buffer = shift;

  # get a buffer-sized chunk of VFs to split and fork on
  my $vfs = $buffer->next();
  return [] unless $vfs && scalar @$vfs;

  my $fork_number = $self->param('fork');
  my $buffer_size = $self->param('buffer_size');
  my $delta = 0.5;
  my $minForkSize = 50;
  my $maxForkSize = int($buffer_size / (2 * $fork_number));
  my $active_forks = 0;
  my (@pids, %by_pid);
  my $sel = IO::Select->new;

  # loop while variants in @$vfs or forks running
  while(@$vfs or $active_forks) {

    # only spawn new forks if we have space
    if($active_forks <= $fork_number) {
      my $numLines = scalar @$vfs;
      my $forkSize = int($numLines / ($fork_number + ($delta * $fork_number)) + $minForkSize ) + 1;

      $forkSize = $maxForkSize if $forkSize > $maxForkSize;

      while(@$vfs && $active_forks <= $fork_number) {

        # create sockets for IPC
        my ($child, $parent);
        socketpair($child, $parent, AF_UNIX, SOCK_STREAM, PF_UNSPEC) or throw("ERROR: Failed to open socketpair: $!");
        $child->autoflush(1);
        $parent->autoflush(1);
        $sel->add($child);

        # readjust forkSize if it's bigger than the remaining buffer
        # otherwise the input buffer will read more from the parser
        $forkSize = scalar @$vfs if $forkSize > scalar @$vfs;
        my @tmp = splice(@$vfs, 0, $forkSize);

        # fork
        my $pid = fork;
        if(!defined($pid)) {
          throw("ERROR: Failed to fork\n");
        }
        elsif($pid) {
          push @pids, $pid;
          $active_forks++;
        }
        elsif($pid == 0) {
          $self->_forked_process($buffer, \@tmp, $parent);
        }
      }
    }

    # read child input
    while(my @ready = $sel->can_read()) {
      my $no_read = 1;

      foreach my $fh(@ready) {
        $no_read++;

        my $line = join('', $fh->getlines());
        next unless $line;
        $no_read = 0;

        my $data = thaw($line);
        next unless $data && $data->{pid};

        # data
        $by_pid{$data->{pid}} = $data->{output} if $data->{output};

        # plugin data
        foreach my $plugin_name(keys %{$data->{plugin_data} || {}}) {
          my ($parent_plugin) = grep {ref($_) eq $plugin_name} @{$self->get_all_Plugins};
          next unless $parent_plugin;

          merge_hashes($parent_plugin, $data->{plugin_data}->{$plugin_name});
        }

        # stats
        merge_hashes($self->stats->{stats}->{counters}, $data->{stats}, 1) if $data->{stats};

        # stderr
        $self->warning_msg($data->{stderr}) if $data->{stderr};

        # finish up
        $sel->remove($fh);
        $fh->close;
        $active_forks--;

        throw("ERROR: Forked process(es) died\n".$data->{die}) if $data->{die};
      }

      # read-through detected, DIE
      throw("\nERROR: Forked process(es) died\n") if $no_read;

      last if $active_forks < $fork_number;
    }
  }

  waitpid($_, 0) for @pids;

  # sort data by dispatched PID order and return
  return [map {@{$by_pid{$_} || []}} @pids];
}

sub _forked_process {
  my $self = shift;
  my $buffer = shift;
  my $vfs = shift;
  my $parent = shift;

  # redirect and capture STDERR
  $self->config->{warning_fh} = *STDERR;
  close STDERR;
  my $stderr;
  open STDERR, '>', \$stderr;

  # reset the input buffer and add a chunk of data to its pre-buffer
  # this way it gets read in on the following next() call
  # which will be made by _buffer_to_output()
  $buffer->{buffer_size} = scalar @$vfs;
  $buffer->reset_buffer();
  $buffer->reset_pre_buffer();
  push @{$buffer->pre_buffer}, @$vfs;

  # reset stats
  $self->stats->{stats}->{counters} = {};

  # reset FASTA DB
  delete($self->config->{_fasta_db});
  $self->fasta_db;

  # we want to capture any deaths and accurately report any errors
  # so we use eval to run the core chunk of the code (_buffer_to_output)
  my $output;
  eval {
    # for testing
    $self->warning_msg('TEST WARNING') if $self->{_test_warning};
    throw('TEST DIE') if $self->{_test_die};

    # the real thing
    $output = $self->_buffer_to_output($buffer);
  };
  my $die = $@;

  # some plugins may cache stuff, check for this and try and
  # reconstitute it into parent's plugin cache
  my $plugin_data;

  foreach my $plugin(@{$self->get_all_Plugins}) {
    next unless $plugin->{has_cache};

    # delete unnecessary stuff and stuff that can't be serialised
    delete $plugin->{$_} for qw(config feature_types variant_feature_types version feature_types_wanted variant_feature_types_wanted params);

    $plugin_data->{ref($plugin)} = $plugin;
  }

  # send everything we've captured to the parent process
  # PID allows parent process to re-sort output to correct order
  print $parent freeze({
    pid => $$,
    output => $output,
    plugin_data => $plugin_data,
    stderr => $stderr,
    die => $die,
    stats => $self->stats->{stats}->{counters},
  });

  exit(0);
}

sub post_setup_checks {
  my $self = shift;

  # disable HGVS if no FASTA file found and it was switched on by --everything
  if(
    $self->param('hgvs') &&
    $self->param('offline') &&
    $self->param('everything') &&
    !$self->fasta_db
  ) {
    $self->status_msg("INFO: Disabling --hgvs; using --offline and no FASTA file found\n");
    $self->param('hgvs', 0);
  }
  
  # offline needs cache, can't use HGVS
  if($self->param('offline')) {
    unless($self->fasta_db) {
      die("ERROR: Cannot generate HGVS coordinates in offline mode without a FASTA file (see --fasta)\n") if $self->param('hgvs');
      die("ERROR: Cannot check reference sequences without a FASTA file (see --fasta)\n") if $self->param('check_ref')
    }
    
    # die("ERROR: Cannot do frequency filtering in offline mode\n") if defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
    die("ERROR: Cannot map to LRGs in offline mode\n") if $self->param('lrg');
  }
    
  # warn user DB will be used for SIFT/PolyPhen/HGVS/frequency/LRG
  if($self->param('cache')) {
        
    # these two def depend on DB
    foreach my $param(grep {$self->param($_)} qw(lrg check_sv)) {
      $self->status_msg("INFO: Database will be accessed when using --$param");
    }

    # and these depend on either DB or FASTA DB
    unless($self->fasta_db) {
      foreach my $param(grep {$self->param($_)} qw(hgvs check_ref)) {
        $self->status_msg("INFO: Database will be accessed when using --$param");
      }
    }
        
    # $self->status_msg("INFO: Database will be accessed when using --check_frequency with population ".$config->{freq_pop}) if defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
  }

  # stats_html should be default, but don't mess if user has already selected one or both
  unless($self->param('stats_html') || $self->param('stats_text')) {
    $self->param('stats_html', 1);
  }

  return 1;
}

sub setup_db_connection {
  my $self = shift;

  return if $self->param('offline');
  return unless $self->param('database') || $self->param('cache');

  # doing this inits the registry and DB connection
  my $reg = $self->registry();

  # check assembly
  if(my $db_assembly = $self->get_database_assembly) {

    my $config_assembly = $self->param('assembly');

    throw(
      "ERROR: Assembly version specified by --assembly (".$config_assembly.
      ") and assembly version in coord_system table (".$db_assembly.") do not match\n".
      (
        $self->param('host') eq 'ensembldb.ensembl.org' ?
        "\nIf using human GRCh37 add \"--port 3337\"".
        " to use the GRCh37 database, or --offline to avoid database connection entirely\n" :
        ''
      )
    ) if $config_assembly && $config_assembly ne $db_assembly;

    # update to database version
    $self->param('assembly', $db_assembly);

    if(!$self->param('assembly')) {
      throw("ERROR: No assembly version specified, use --assembly [version] or check the coord_system table in your core database\n");
    }
  }

  # update species, e.g. if user has input "human" we get "homo_sapiens"
  $self->species($reg->get_alias($self->param('species')));

  return 1;
}

sub get_all_AnnotationSources {
  my $self = shift;

  if(!exists($self->{_annotation_sources})) {
    my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $self->config});
    $self->{_annotation_sources} = $asa->get_all;
  }

  return $self->{_annotation_sources};
}

sub get_Parser {
  my $self = shift;

  if(!exists($self->{parser})) {

    # user given input data as string (REST)?
    if(my $input_data = $self->param('input_data')) {
      open IN, '<', \$input_data;
      $self->param('input_file', *IN);
    }

    $self->{parser} = Bio::EnsEMBL::VEP::Parser->new({
      config            => $self->config,
      format            => $self->param('format'),
      file              => $self->param('input_file'),
      valid_chromosomes => $self->get_valid_chromosomes,
    })
  }

  return $self->{parser};
}

sub get_InputBuffer {
  my $self = shift;

  if(!exists($self->{input_buffer})) {
    $self->{input_buffer} = Bio::EnsEMBL::VEP::InputBuffer->new({
      config => $self->config,
      parser => $self->get_Parser
    });
  }

  return $self->{input_buffer};
}

sub get_OutputFactory {
  my $self = shift;

  if(!exists($self->{output_factory})) {
    $self->{output_factory} = Bio::EnsEMBL::VEP::OutputFactory->new({
      config      => $self->config,
      format      => $self->param('output_format'),
      header_info => $self->get_output_header_info,
      plugins     => $self->get_all_Plugins,
    });
  }

  return $self->{output_factory};
}

sub get_all_Plugins {
  my $self = shift;

  if(!exists($self->{plugins})) {
    my @plugins = ();

    unshift @INC, $self->param('dir_plugins') || $self->param('dir').'/Plugins';

    PLUGIN: foreach my $plugin_config(@{$self->param('plugin') || []}) {

      # parse out the module name and parameters
      my ($module, @params) = split /,/, $plugin_config;

      # check we can use the module      
      eval qq{
        use $module;
      };
      if($@) {
        my $msg = "Failed to compile plugin $module: $@\n";
        throw($msg) if $self->param('safe');
        $self->warning_msg($msg);
        next;
      }
      
      # now check we can instantiate it, passing any parameters to the constructor      
      my $instance;
      
      eval {
        $instance = $module->new($self->config, @params);
      };
      if($@) {
        my $msg = "Failed to instantiate plugin $module: $@\n";
        throw($msg) if $self->param('safe');
        $self->warning_msg($msg);
        next;
      }

      # check that the versions match
      
      #my $plugin_version;
      #
      #if ($instance->can('version')) {
      #  $plugin_version = $instance->version;
      #}
      #
      #my $version_ok = 1;
      #
      #if ($plugin_version) {
      #  my ($plugin_major, $plugin_minor, $plugin_maintenance) = split /\./, $plugin_version;
      #  my ($major, $minor, $maintenance) = split /\./, $VERSION;
      #
      #  if ($plugin_major != $major) {
      #    debug("Warning: plugin $plugin version ($plugin_version) does not match the current VEP version ($VERSION)") unless defined($config->{quiet});
      #    $version_ok = 0;
      #  }
      #}
      #else {
      #  debug("Warning: plugin $plugin does not define a version number") unless defined($config->{quiet});
      #  $version_ok = 0;
      #}
      #
      #debug("You may experience unexpected behaviour with this plugin") unless defined($config->{quiet}) || $version_ok;

      # check that it implements all necessary methods
      
      for my $required(qw(run get_header_info check_feature_type check_variant_feature_type feature_types)) {
        unless($instance->can($required)) {
          my $msg = "Plugin $module doesn't implement a required method '$required', does it inherit from BaseVepPlugin?\n";
          throw($msg) if $self->param('safe');
          $self->warning_msg($msg);
          next PLUGIN;
        }
      }
       
      # all's good, so save the instance in our list of plugins      
      push @plugins, $instance;
      
      $self->status_msg("Loaded plugin: $module");

      # for convenience, check if the plugin wants regulatory stuff and turn on the config option if so
      if (grep { $_ =~ /motif|regulatory/i } @{ $instance->feature_types }) {
        $self->status_msg("Fetching regulatory features for plugin: $module");
        $self->param('regulatory', 1);
      }
    }

    $self->{plugins} = \@plugins;
  }

  return $self->{plugins};
}

sub get_output_file_handle {
  my $self = shift;

  unless(exists($self->{output_file_handle})) {

    my $output_file_handle = FileHandle->new();

    my $output_file_name = $self->param('output_file');
      
    # check if file exists
    if(-e $output_file_name && !$self->param('force_overwrite')) {
      throw("ERROR: Output file $output_file_name already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }
    
    if(uc($output_file_name) eq 'STDOUT') {
      $output_file_handle = *STDOUT;
    }
    else {
      $output_file_handle->open(">$output_file_name") or throw("ERROR: Could not write to output file $output_file_name\n");
    }

    $self->{output_file_handle} = $output_file_handle;
  }

  return $self->{output_file_handle};
}

sub get_stats_file_handle {
  my $self = shift;
  my $type = shift;

  my $stats_file_root = $self->param('stats_file') || $self->param('output_file').'_summary';

  my $file_name;

  if($stats_file_root =~ m/^(.+?)\.(.+)$/) {
    my ($root, $ext) = ($1, $2);
    if(lc($ext) eq $type) {
      $file_name = $stats_file_root;
    }
    elsif(lc($ext) =~ /^(txt|html)$/) {
      $file_name = $root.'.'.$type;
    }
    else {
      $file_name = $stats_file_root.'.'.$type;
    }
  }
  else {
    $file_name = $stats_file_root.'.'.$type;
  }
      
  # check if file exists
  if(-e $file_name && !$self->param('force_overwrite')) {
    throw("ERROR: Stats file $file_name already exists. Specify a different output file with --stats_file or overwrite existing file with --force_overwrite\n");
  }
  
  my $fh = FileHandle->new();
  $fh->open(">$file_name") or throw("ERROR: Could not write to stats file $file_name\n");

  return $fh;
}

sub dump_stats {
  my $self = shift;

  unless($self->param('no_stats')) {

    if($self->param('stats_text')) {
      my $fh = $self->get_stats_file_handle('txt');
      $self->stats->dump_text($fh);
      close $fh;
    }

    if($self->param('stats_html')) {
      my $fh = $self->get_stats_file_handle('html');
      $self->stats->dump_html($fh);
      close $fh;
    }
  }
}

sub get_output_header_info {
  my $self = shift;

  if(!exists($self->{output_header_info})) {

    my $info = {
      time          => get_time,
      vep_version   => $Bio::EnsEMBL::VEP::Constants::VERSION,
      api_version   => $self->registry->software_version,
      input_headers => $self->get_Parser->headers,
    };

    if(my $mca = $self->get_adaptor('core', 'MetaContainer')) {
      $info->{db_name} = $mca->dbc->dbname;
      $info->{db_host} = $mca->dbc->host;
      $info->{db_version} = $mca->get_schema_version;
    }

    foreach my $as(@{$self->get_all_AnnotationSources}) {
      my $as_info = $as->info;
      $info->{version_data}->{$_} ||= $as_info->{$_} for keys %$as_info;
      $info->{cache_dir} ||= $as->dir if $as->can('dir');
      push @{$info->{custom_info}}, $as_info->{custom_info} if $as_info->{custom_info};
    }

    $self->{output_header_info} = $info;
  }

  return $self->{output_header_info};
}

sub get_valid_chromosomes {
  my $self = shift;

  if(!exists($self->{valid_chromosomes})) {

    my %valid = ();

    foreach my $as(@{$self->get_all_AnnotationSources}) {
      next unless $as->can('get_valid_chromosomes');
      $valid{$_} = 1 for @{$as->get_valid_chromosomes};
    }

    $self->{valid_chromosomes} = [sort keys %valid];
  }

  return $self->{valid_chromosomes};
}

1;