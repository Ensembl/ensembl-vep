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

# EnsEMBL module for Bio::EnsEMBL::VEP::Runner
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Runner - runner class for VEP

=head1 SYNOPSIS

my $runner = Bio::EnsEMBL::VEP::Runner->new({
  input_file    => 'variants.vcf',
  output_file   => 'variants_annotated.vcf',
  output_format => 'vcf',
  cache         => 1,
});

$runner->param('fork', 4);

$runner->run();

=head1 DESCRIPTION

The Runner class is the entry point for access to most
of VEP's functionality.

Creating and configuring Runner should be all most users
need to access VEP functionality.

Shares a base class, BaseRunner, with Haplo::Runner.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Runner;

use base qw(Bio::EnsEMBL::VEP::BaseRunner);

use Storable qw(freeze thaw);
use IO::Socket;
use IO::Select;
use FileHandle;
use Scalar::Util qw(openhandle);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(get_time merge_hashes);
use Bio::EnsEMBL::VEP::Constants;
use Bio::EnsEMBL::VEP::Parser;
use Bio::EnsEMBL::VEP::InputBuffer;
use Bio::EnsEMBL::VEP::OutputFactory;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(revert_fasta);


=head2 init

  Example    : $runner->init()
  Description: Runs initialisation and setup for the runner. This includes:
                - initialising stats
                - setting up database connection
                - initialising plugins
                - getting all AnnotationSources
                - loading chromosome synonyms
                - setting up FASTA file for sequence access
                - getting InputBuffer object
                - running post-setup checks
               Runs only once - subsequent calls are not executed.
  Returntype : bool
  Exceptions : none
  Caller     : run(), run_rest(), next_output_line()
  Status     : Stable

=cut

sub init {
  my $self = shift;

  return 1 if $self->{_initialized};

  # log start time
  $self->stats->start_time();

  if(defined($self->param('var_synonyms')) && !$self->param('check_existing')) {
    $self->status_msg("INFO: Colocated variants will be provided when using --var_synonyms\n");
    $self->param('check_existing', 1);
  }

  # setup DB connection
  $self->setup_db_connection();

  my $plugins = $self->get_all_Plugins();

  # get all annotation sources
  my $annotation_sources = $self->get_all_AnnotationSources();

  # get chromosome synoyms
  $self->chromosome_synonyms($self->param('synonyms'));

  # setup FASTA file DB
  $self->fasta_db();

  my $buffer = $self->get_InputBuffer();

  $self->post_setup_checks();

  $self->stats->info($self->get_output_header_info);

  # get stats file handle now
  # otherwise we don't know if this will fail til the end of the run
  unless($self->param('no_stats')) {
    $self->stats->log_chromosomes($self->chr_lengths);
    $self->get_stats_file_handle('txt') if $self->param('stats_text');
    $self->get_stats_file_handle('html') if $self->param('stats_html');
  }

  return $self->{_initialized} = 1;
}


=head2 finish

  Example    : $runner->finish()
  Description: Finalises run and restores some states. This includes:
                - calling plugins' finish() method
                - dumping VEP run statistics
                - resetting FASTA file config
  Returntype : none
  Exceptions : none
  Caller     : run(), run_rest()
  Status     : Stable

=cut

sub finish {
  my $self = shift;

  foreach my $plugin(@{$self->get_all_Plugins}) {
    $plugin->finish() if $plugin->can('finish');
  }

  $self->dump_stats unless $self->param('no_stats');

  revert_fasta() if $self->fasta_db;

  if(my $warning_fh = $self->config->{warning_fh}) {
    close($warning_fh) if openhandle($warning_fh);
  }
}


=head2 run

  Example    : $runner->run()
  Description: Run VEP: initialise, get and write output, finish.
  Returntype : bool
  Exceptions : none
  Caller     : vep
  Status     : Stable

=cut

sub run {
  my $self = shift;

  $self->init();

  my $fh = $self->get_output_file_handle();

  unless($self->param('no_headers')) {
    print $fh "$_\n" for @{$self->get_OutputFactory->headers};
  }

  while(my $line = $self->next_output_line) {
    print $fh "$line\n";
  }

  close $fh;

  $self->finish();

  return 1;
}


=head2 run_rest

  Arg 1      : (optional) string $input
  Example    : $results = $runner->run_rest($input)
  Description: Run VEP as used from REST API: initialise, get output, finish,
               return output as listref of hashrefs.
  Returntype : listref of hashrefs
  Exceptions : none
  Caller     : REST API VEP Controller
  Status     : Stable

=cut

sub run_rest {
  my $self = shift;
  my $input = shift;

  $self->internalise_warnings();

  $self->param('input_data', $input);
  $self->param('output_format', 'json');
  $self->param($_, 1) for qw(no_escape safe quiet no_stats);

  $self->init();

  my @return = ();

  while(my $hash = $self->next_output_line(1)) {
    push @return, $hash;
  }

  $self->finish();

  close WARNINGS;

  return \@return;
}


=head2 internalise_warnings

  Example    : $runner->internalise_warnings();
               $runner->run($dodgy_input);
               $warnings = $runner->warnings();
  Description: Instead of printing warnings/errors to STDERR (or file
               if using --warning_file), store internally. Warnings
               can then be retrieved using warnings();
  Returntype : none
  Exceptions : none
  Caller     : run_rest()
  Status     : Stable

=cut

sub internalise_warnings {
  my $self = shift;
  $self->{_warning_string} = '';
  open WARNINGS, '>', \$self->{_warning_string};
  $self->config->{warning_fh} = *WARNINGS;
}


=head2 warnings

  Example    : $results = $runner->run_rest($dodgy_input);
               $warnings = $runner->warnings();
  Description: Used after run_rest() to check for warnings/errors
               returns an arrayref of hashrefs with the following structure:
               {
                 type => 'ERROR',
                 msg => 'Something exploded',
                 stack => 'STACK Foo::Bar::foobar()\nSTACK Boo::Far::moocar()'
               }
  Returntype : listref of hashrefs
  Exceptions : none
  Caller     : REST API VEP Controller
  Status     : Stable

=cut

sub warnings {
  my $self = shift;

  my @warnings;

  my %current_warning = ();

  if(my $warnings = $self->{_warning_string}) {
    foreach my $line(split /\n+/, $warnings) {

      if($line =~ /^(WARNING|ERROR)\s*\:\s*(.+)$/) {

        if(keys %current_warning) {
          my %copy = %current_warning;
          push @warnings, \%copy;
        }

        %current_warning = (
          type => $1,
          msg => $2,
        );
      }
      elsif($line =~ /^MSG\s*\:\s*(.+)$/) {
        $current_warning{msg} .= ': '.$1;
      }
      elsif($line =~ /^STACK/) {
        $current_warning{stack} .= $line."\n";
      }
    }
  }

  if(keys %current_warning) {
    my %copy = %current_warning;
    push @warnings, \%copy;
  }

  return \@warnings;
}


=head2 next_output_line

  Arg 1      : (optional) bool $output_as_hash
  Example    : $line = $runner->next_output_line();
  Description: Get the next line of output. If $output_as_hash is set
               to a true value, a hashref of data is returned; otherwise
               return value is a string.
  Returntype : string or hashref
  Exceptions : none
  Caller     : run()
  Status     : Stable

=cut

sub next_output_line {
  my $self = shift;
  my $output_as_hash = shift;

  my $output_buffer = $self->{_output_buffer} ||= [];

  return shift @$output_buffer if @$output_buffer;

  $self->init();

  $self->_set_package_variables();

  if($self->param('fork')) {
    push @$output_buffer, @{$self->_forked_buffer_to_output($self->get_InputBuffer, $output_as_hash)};
  }
  else {
    push @$output_buffer, @{$self->_buffer_to_output($self->get_InputBuffer, $output_as_hash)};
  }

  $self->_reset_package_variables();

  return @$output_buffer ? shift @$output_buffer : undef;
}


=head2 _buffer_to_output

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Arg 2      : (optional) bool $output_as_hash
  Example    : $output = $runner->_buffer_to_output($ib);
  Description: Get the next buffer "fill" and run annotate. This may be
               delegated to a number of forked processes or run in this
               process depending on configuration. If $output_as_hash is set
               to a true value, arrayrefs of hashref of data are returned; otherwise
               return value is an arrayref of strings.
  Returntype : arrayref of strings or hashrefs
  Exceptions : none
  Caller     : next_output_line(), _forked_process()
  Status     : Stable

=cut

sub _buffer_to_output {
  my $self = shift;
  my $input_buffer = shift;
  my $output_as_hash = shift;

  my @output;
  my $vfs = $input_buffer->next();
  return [] unless $vfs && scalar @$vfs;

  my $output_factory = $self->get_OutputFactory;

  foreach my $as(@{$self->get_all_AnnotationSources}) {
    $as->annotate_InputBuffer($input_buffer);
  }
    
  $input_buffer->finish_annotation;

  # if the input buffer got reduced to nothing we have to re-run this method here
  # this can happen if an annotation source filters out variants
  # NB we dont want to do this if we've been called from within a fork otherwise duplication will happen
  return $self->_buffer_to_output($input_buffer, $output_as_hash) unless scalar @{$input_buffer->buffer} || $self->param('fork');

  if($output_as_hash) {
    push @output, @{$output_factory->get_all_output_hashes_by_InputBuffer($input_buffer)};
  }
  else {
    push @output, @{$output_factory->get_all_lines_by_InputBuffer($input_buffer)};
  }
  # if output is empty here rerun because there could still be VFs left coming from the parser 
  return $self->_buffer_to_output($input_buffer, $output_as_hash) unless scalar @output  || $self->param('fork');

  return \@output;
}


=head2 _forked_buffer_to_output

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Arg 2      : (optional) bool $output_as_hash
  Example    : $output = $runner->_forked_buffer_to_output($ib);
  Description: Fork-based wrapper for _buffer_to_output. Divides variants
               from InputBuffer into chunks and forks off a process to
               run each chunk. Output is collected from socketpairs that
               allow communication between the forked children and this
               process; data is written to sockets as hashes serialised
               with Storable.
  Returntype : arrayref of strings or hashrefs
  Exceptions : none
  Caller     : next_output_line()
  Status     : Stable

=cut

sub _forked_buffer_to_output {
  my $self = shift;
  my $buffer = shift;
  my $output_as_hash = shift;

  # get a buffer-sized chunk of VFs to split and fork on
  my $vfs = $buffer->next();
  return [] unless $vfs && scalar @$vfs;

  # log lines read here
  $self->stats->log_lines_read($buffer->parser->line_number) if $buffer->parser && !$self->param('no_stats');

  my $fork_number = $self->param('fork');
  my $buffer_size = $self->param('buffer_size');
  my $delta = 0.5;
  my $minForkSize = 50;
  my $maxForkSize = int($buffer_size / (2 * $fork_number)) || 1;
  my $active_forks = 0;
  my (@pids, %by_pid, %children_reported_back);
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

          # for testing
          kill $self->{_kill_child}, $pid if $self->{_kill_child};
        }
        elsif($pid == 0) {
          $self->_forked_process($buffer, \@tmp, $parent, $output_as_hash);
        }
      }
    }

    # read child output
    while(my @ready = $sel->can_read()) {
      my $no_read = 1;

      foreach my $fh(@ready) {
        $no_read++;

        my $line = join('', $fh->getlines());
        next unless $line;
        $no_read = 0;

        my $data = thaw($line);
        next unless $data && $data->{pid};

        # forked process died
        die(sprintf("%sDied in forked process %i\n", $data->{die}, $data->{pid})) if $data->{die};

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
        $self->warning_msg($data->{pid}." : ".$data->{stderr}) if $data->{stderr};

        # merge in cached data
        if(my $fork_vep_cache = $data->{_VEP_CACHE}) {
          my $cache = $main::_VEP_CACHE ||= {};
          merge_hashes($cache, $fork_vep_cache);
        }

        # report this child as done
        $children_reported_back{$data->{pid}} = 1;

        # finish up
        $sel->remove($fh);
        $fh->close;
        $active_forks--;
      }

      # read-through detected, DIE
      throw("\nERROR: Forked process(es) died: read-through of cross-process communication detected\n") if $no_read;

      last if $active_forks < $fork_number;
    }
  }

  for my $pid(@pids) {
    if(waitpid($pid, 0) > 0) {
      my ($rc, $sig, $core) = ($? >> 8, $? & 127, $? & 128);
      if($core){
        throw("\nERROR: Forked process PID $pid dumped core\n");
      }
      elsif($sig == 9) {
        throw("\nERROR: Forked process PID $pid was killed\n");
      }
      elsif($rc != 0) {
        throw("\nERROR: Forked process PID $pid returned non-zero return code $rc\n");
      }
    }
    else {
      throw("\nERROR: Forked process PID $pid disappeared\n");
    }
  }

  # check any children didn't report back
  my @unclean = grep {!$children_reported_back{$_}} @pids;
  throw("\nERROR: Forked process(es) ".join(", ", @unclean)." did not finish cleanly") if @unclean;

  # sort data by dispatched PID order
  my @return = map {@{$by_pid{$_} || []}} @pids;

  # if the input buffer got reduced to nothing we have to re-run this method here
  # this can happen if an annotation source filters out variants
  return $self->_forked_buffer_to_output($buffer, $output_as_hash) unless scalar @return;

  return \@return;
}


=head2 _forked_process

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Arg 2      : arrayref of Bio::EnsEMBL::VariationFeature $vfs
  Arg 3      : filehandle $parent_socket
  Arg 4      : bool $output_as_hash
  Example    : $runner->_forked_process($ib, $vfs, $socket);
  Description: Forked process is run within processes forked off by
               _forked_buffer_to_output(). It essentially executes
               _buffer_to_output() on the chunk of input given, and
               prints output as hashrefs serialised with Storable to
               a socketpair that is read by the parent process. These
               serialised hashrefs also contain stats, cached plugin
               data and any warnings or die states.
  Returntype : none (exits process)
  Exceptions : none
  Caller     : _forked_buffer_to_output()
  Status     : Stable

=cut

sub _forked_process {
  my $self = shift;
  my $buffer = shift;
  my $vfs = shift;
  my $parent = shift;
  my $output_as_hash = shift;

  # for testing
  kill $self->{_kill_self}, $$ if $self->{_kill_self};

  # simulate a memory leak
  # my $loop = 1;
  # my $mem_leak = sub { my ($a, $b); $a = \$b; $b = \$a; print STDERR "$$ ".($loop / 1e6)."\n" if ++$loop % 1e6 == 0};
  # &$mem_leak() while(1);

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

  # reset custom sources' parsers
  # otherwise we get cross-pollution between forks reading from the same filehandles (I think)
  delete $_->{parser} for @{$self->get_all_AnnotationSources};

  # we want to capture any deaths and accurately report any errors
  # so we use eval to run the core chunk of the code (_buffer_to_output)
  my $output;
  eval {
    # for testing
    $self->warning_msg('TEST WARNING') if $self->{_test_warning};
    throw('TEST DIE') if $self->{_test_die};

    # the real thing
    $output = $self->_buffer_to_output($buffer, $output_as_hash);
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
    _VEP_CACHE => $main::_VEP_CACHE,
  });

  exit(0);
}


=head2 post_setup_checks

  Example    : $ok = $runner->post_setup_checks();
  Description: Run a set of checks after initialisation:
                - disable HGVS if using --offline with no FASTA
                - check for set options incompatible with --offline
                - warn user certain flags will access DB with --cache
  Returntype : bool
  Exceptions : throws if:
                - --hgvs and --offline enabled without a FASTA file
                - --check_ref and --offline enabled without a FASTA file
                - --lrg and --offline set
  Caller     : init()
  Status     : Stable

=cut

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
    $self->param($_, 0) for qw(hgvs hgvsc hgvsp);
  }
  
  # offline needs cache, can't use HGVS
  if($self->param('offline')) {
    unless($self->fasta_db) {
      throw("ERROR: Cannot generate HGVS coordinates (--hgvs) in offline mode without a FASTA file (see --fasta)\n") if $self->param('hgvs') || $self->param('hgvsc') || $self->param('hgvsp');
      throw("ERROR: Cannot check reference sequences (--check_ref) without a FASTA file (see --fasta)\n") if $self->param('check_ref');
      throw("ERROR: Cannot look up reference sequences (--lookup_ref) without a FASTA file (see --fasta)\n") if $self->param('lookup_ref');
      throw("ERROR: Cannot use transcript reference sequences (--use_transcript_ref) without a FASTA file (see --fasta); you may wish to use --use_given_ref\n") if $self->param('use_transcript_ref');
    }
    
    # throw("ERROR: Cannot do frequency filtering in offline mode\n") if defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
    throw("ERROR: Cannot map to LRGs in offline mode\n") if $self->param('lrg');
  }
    
  # warn user DB will be used for SIFT/PolyPhen/HGVS/frequency/LRG
  if($self->param('cache')) {
        
    # these two def depend on DB
    foreach my $param(grep {$self->param($_)} qw(lrg check_sv)) {
      $self->status_msg("INFO: Database will be accessed when using --$param");
    }

    # and these depend on either DB or FASTA DB
    unless($self->fasta_db) {
      foreach my $param(grep {$self->param($_)} qw(hgvs hgvsc hgvsp lookup_ref check_ref use_transcript_ref)) {
        $self->status_msg("INFO: Database will be accessed when using --$param");
      }
    }
        
    # $self->status_msg("INFO: Database will be accessed when using --check_frequency with population ".$config->{freq_pop}) if defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
  }

  # stats_html should be default, but don't mess if user has already selected one or both
  unless($self->param('stats_html') || $self->param('stats_text')) {
    $self->param('stats_html', 1);
  }
  
  if(defined($self->param('shift_hgvs')) && $self->param('shift_hgvs') eq 0) {
    $self->status_msg("INFO: --shift_hgvs has been set to zero, setting shift_3prime to 0\n");
    $self->param('shift_3prime', 0);
  }

  return 1;
}


=head2 get_Parser

  Example    : my $parser = $runner->get_Parser();
  Description: Gets the Parser configured for this runner.
  Returntype : Bio::EnsEMBL::VEP::Parser
  Exceptions : none
  Caller     : getInputBuffer()
  Status     : Stable

=cut

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
      valid_chromosomes => $self->valid_chromosomes,
    })
  }

  return $self->{parser};
}


=head2 get_InputBuffer

  Example    : my $ib = $runner->get_InputBuffer();
  Description: Gets the InputBuffer configured for this runner.
  Returntype : Bio::EnsEMBL::VEP::InputBuffer
  Exceptions : none
  Caller     : init(), run()
  Status     : Stable

=cut

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


=head2 get_OutputFactory

  Example    : my $ib = $runner->get_OutputFactory();
  Description: Gets the OutputFactory configured for this runner.
  Returntype : Bio::EnsEMBL::VEP::OutputFactory
  Exceptions : none
  Caller     : init(), run()
  Status     : Stable

=cut

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


=head2 get_all_Plugins

  Example    : my $plugins = $runner->get_all_Plugins();
  Description: Takes user config for plugins and initialises plugin object
               instances. Failure to instantiate a plugin instance will not
               throw by default; this must be enabled with --safe (as used
               by web VEP and REST API).
  Returntype : arrayref of Bio::EnsEMBL::Variation::Utils::BaseVepPlugin
  Exceptions : none
  Caller     : init()
  Status     : Stable

=cut

sub get_all_Plugins {
  my $self = shift;

  if(!defined($self->{plugins})) {
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
        $instance = $module->new($self->config->{_params}, @params);
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
      
      # $self->status_msg("Loaded plugin: $module");

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

1;
