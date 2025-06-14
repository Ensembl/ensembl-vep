=head1 LICENSE

Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Base class used by all Bio::EnsEMBL::VEP::* classes

Contains methods for accessing common functionality.

=head1 METHODS

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
use File::Spec;
use FileHandle;

use base qw(Exporter);

our @EXPORT_OK = qw(_fetch_chr_synonyms);

=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config,
                 arg1   => val1,
                 ... 
               }
  Example    : $obj = Bio::EnsEMBL::VEP::*->new($args);
  Description: Constructor method used by most VEP classes
  Returntype : Bio::EnsEMBL::VEP::BaseVEP
  Exceptions : throws if config arg is not Bio::EnsEMBL::VEP::Config
  Caller     : general
  Status     : Stable

=cut

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


=head2 config

  Example    : $config = $obj->config()
  Description: Gets associated Bio::EnsEMBL::VEP::Config
  Returntype : Bio::EnsEMBL::VEP::Config
  Exceptions : none
  Caller     : param()
  Status     : Stable

=cut

sub config {
  return $_[0]->{_config};
}


=head2 param

  Examples   : $value = $obj->param($key)
               $obj->param($key, $value)
  Description: Getter/setter for config parameters
  Returntype : string (typically)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub param {
  my $self = shift;
  return $self->config->param(@_);
}


=head2 stats

  Example    : $stats = $obj->stats()
  Description: Get Bio::EnsEMBL::VEP::Stats object for this VEP run. Cached
               on shared config object so all derived objects will share.
  Returntype : Bio::EnsEMBL::VEP::Stats
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stats {
  my $self = shift;

  my $config = $self->config;

  if(!exists($config->{_stats})) {
    $config->{_stats} = Bio::EnsEMBL::VEP::Stats->new({config => $self->config});
  }

  return $config->{_stats};
}


=head2 skipped_variants

  Arg 1      : (optional) hash $value with line, line_number and description
  Example    : $skipped_variants = $obj->skipped_variants()
  Description: Get/set skipped variants for this VEP run. Cached on shared
               config object so all derived objects will share.
  Returntype : Arrayref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub skipped_variants {
  my $self = shift;
  my $value = shift;

  my $config = $self->config;

  if (!exists($config->{_skipped_variants})) {
    $config->{_skipped_variants} = [];
  }

  if (defined $value) {
    push @{ $config->{_skipped_variants} }, $value;
  }

  return $config->{_skipped_variants};
}


=head2 species

  Example    : $species = $obj->species()
  Description: Getter/setter species name configured for this VEP run
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

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


=head2 registry

  Arg 1      : (optional) Bio::EnsEMBL::Registry $registry
  Example    : $registry = $obj->registry()
  Description: Invokes database connection via Bio::EnsEMBL::Registry.
               Can also be set as first arg; this is used by REST.
  Returntype : Bio::EnsEMBL::Registry
  Exceptions : none
  Caller     : setup_db_connection(), REST
  Status     : Stable

=cut

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
        $registry_file = File::Spec->rel2abs($registry_file);
        throw("ERROR: Registry file $registry_file not found") unless -e $registry_file;
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

        my $host       = $self->param('host');
        my $user       = $self->param('user');
        my $port       = $self->param('port');
        my $db_version = $self->param('db_version');
        $reg->load_registry_from_db(
          -host       => $host,
          -user       => $user,
          -pass       => $self->param('password'),
          -port       => $self->param('port'),
          -db_version => $db_version,
          -species    => $species,
          -verbose    => $self->param('verbose'),
          -no_cache   => $self->param('no_slice_cache'),
        );

        my @db_names = map { $_->dbc->dbname } @{ $reg->get_all_DBAdaptors };
        my $hostname = ($user eq 'anonymous' ? '' : $user . '@') . "$host:$port";
        $self->warning_msg("No database names in $hostname contain version $db_version")
          if defined $db_version and not grep { /$db_version/ } @db_names;
      }

      eval { $reg->set_reconnect_when_lost() };

      # silence version check warnings, we'll frequently use non-current DB versions
      $reg->no_version_check(1);
    }

    $config->{_registry} = $reg;
  }

  # copy to raw_config for plugin backward compatibility
  $config->{_params}->{reg} ||= $config->{_registry};

  return $config->{_registry};
}


=head2 get_adaptor

  Arg 1      : string $group (core, variation, funcgen)
  Arg 2      : string $adaptor_type
  Example    : $slice_adaptor = $obj->get_adaptor('core', 'slice')
  Description: Gets database adaptors for given group and type. Uses internal
               cache for speedup, bypasses requirement for giving species.
               Also fetches "fake" adaptors with no DB connection where necessary.
  Returntype : Bio::EnsEMBL::BaseAdaptor
  Exceptions : throws if no group and/or type given
  Caller     : general
  Status     : Stable

=cut

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
      $ad = $self->registry->get_adaptor($self->species, $group, $type);
    }

    $ad ||= $self->_get_fake_adaptor($group, $type);

    $cache->{$group}->{$type} = $ad;
  }

  return $cache->{$group}->{$type};
}


=head2 _get_fake_adaptor

  Arg 1      : string $group (core, variation, funcgen)
  Arg 2      : string $adaptor_type
  Example    : $fake_vf_adaptor = $obj->get_adaptor('variation', 'VariationFeature')
  Description: Fetches "fake" adaptors with no DB connection.
  Returntype : Bio::EnsEMBL::BaseAdaptor
  Exceptions : none
  Caller     : get_adaptor()
  Status     : Stable

=cut

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


=head2 get_slice

  Arg 1      : string $chr
  Example    : $slice = $obj->get_slice('MT')
  Description: Gets whole-chromosome slices from slice adaptor. Creates fake slices
               using FASTA database if no database connection available.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_slice {
  my $self = shift;
  my $orig_chr = shift;

  my $chr = $self->get_source_chr_name($orig_chr, 'slices', [keys %{$self->chr_lengths}]);

  my $cache = $self->{_slice_cache} ||= {};

  if(!exists($cache->{$chr})) {
    my $slice;

    if(my $sa = $self->get_adaptor(($self->{core_type} || $self->param('core_type')), 'Slice')) {
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


=head2 get_database_assembly

  Example    : $assembly = $obj->get_database_assembly()
  Description: Get the assembly version for the connected core database.
  Returntype : string
  Exceptions : none
  Caller     : BaseRunner, AnnotationSource::Database::Transcript
  Status     : Stable

=cut

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


=head2 fasta_db

  Example    : $fasta_db = $obj->fasta_db()
  Description: Sets up FASTA file handler using Bio::EnsEMBL::Variation::Utils::FastaSequence
  Returntype : Bio::DB::HTS::Faidx or Bio::DB::Fasta
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fasta_db {
  my $self = shift;

  if(!exists($self->config->{_fasta_db})) {
    my $fasta_db;

    # work out fasta path if fasta_dir given
    if(my $fasta_dir = $self->param('fasta_dir')) {
      if(opendir DIR, $fasta_dir) {
        my ($species, $assembly) = ($self->species, $self->param('assembly'));
        my $can_use_faidx = $Bio::EnsEMBL::Variation::Utils::FastaSequence::CAN_USE_FAIDX;

        my @files =
          sort {  # prefer toplevel, bgzipped
            $a !~ /toplevel/ <=> $b !~ /toplevel/ ||
            $a !~ /gz$/ <=> $b !~ /gz$/
          }
          grep {  # find only those matching species, assembly, bgzipped only if Bio::DB::HTS::Faidx installed
            (/\.fa$/ || ($can_use_faidx && /\.fa.gz$/)) &&
            /$species/i &&
            /$assembly/i
          } readdir DIR;
        closedir DIR;

        if(@files) {
          $self->param('fasta', $fasta_dir.'/'.$files[0]);
        }
        else {
          $self->warning_msg("No compatible FASTA files found in $fasta_dir");
        }
      }
    }

    if(my $fasta_file = $self->param('fasta')) {

      $fasta_db = setup_fasta(
        -FASTA => $fasta_file,
        -ASSEMBLY => $self->param('assembly'),
        -OFFLINE => $self->param('offline'),
        -SYNONYMS => $self->chromosome_synonyms,
      );
    }

    $self->config->{_fasta_db} = $fasta_db;
  }

  return $self->config->{_fasta_db};
}


=head2 chr_lengths

  Example    : $chr_lengths = $obj->chr_lengths()
  Description: Gets all valid chromosome names and their lengths from slice
               adaptor if available, or FASTA database if not
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub chr_lengths {
  my $self = shift;

  if(!exists($self->{_chr_lengths})) {
    my %chr_lengths;
    
    if(
      my $sa = $self->get_adaptor(
        $self->{core_type} || $self->param('core_type'),
        'Slice'
      )
    ) {
      foreach my $slice(@{$sa->fetch_all('toplevel')}, @{$sa->fetch_all('lrg', undef, 1, undef, 1)}) {
        my $chr = $slice->seq_region_name;
        $chr_lengths{$chr} = $slice->length;
      }
    }

    elsif(my $fasta_db = $self->fasta_db) {
      %chr_lengths =
        map {$_ => $fasta_db->length($_)}
        $fasta_db->isa('Bio::DB::Fasta') ?
        (grep {$_ !~ /^\_\_/} $fasta_db->get_all_primary_ids) :
        ($fasta_db->get_all_sequence_ids);
    }

    $self->{_chr_lengths} = \%chr_lengths;
  }

  return $self->{_chr_lengths};
}


=head2 _fetch_chr_synonyms

  Arg 1      : string $SeqRegionSynonyms_adaptor
  Arg 2      : string $slice_adaptor
  Example    : $tree = $obj->_fetch_chr_synonyms($srsa, $sa)
  Description: Fetches chromosome synonyms
  Returntype : hashref
  Exceptions : none
  Caller     : chromosome_synonyms()
  Status     : Stable

=cut

sub _fetch_chr_synonyms {
  my ($srsa, $sa) = @_;

  # synonyms can be indirect i.e. A <-> B <-> C
  # and there may not be a direct link between A <-> C in the DB
  # so let's allow for one level of indirection
  my $tree = {};

  my @all_syns = @{$srsa->fetch_all};

  # To prevent memory errors, species with many seq_regions will not search for indirect synonyms
  # With 4GB memory supplied, errors seem to start with species with ~50k seq regions, so 
  # we've chosen 40k as the limit for calculating indirect synonyms 
  my $syn_threshold = 40000; 
  foreach my $syn(@all_syns) {
    my $syn_slice = $sa->fetch_by_seq_region_id($syn->seq_region_id);
    next unless $syn_slice;
    my ($a, $b) = sort ($syn_slice->seq_region_name, $syn->name);
    $tree->{$a}->{$b} = 1;
    unless(scalar(@all_syns) > $syn_threshold){ 
      $tree->{$_}->{$b} = 1 for keys %{$tree->{$a} || {}};
      $tree->{$_}->{$a} = 1 for keys %{$tree->{$b} || {}};
    }
  }
  return $tree;
}


=head2 chromosome_synonyms
  
  Arg 1      : (optional) string $synonyms_file
  Example    : $syns = $obj->chromosome_synonyms()
  Description: Gets hashref of chromosome synonyms; synonyms from $synonyms_file
               are added if specfied. 
  Returntype : hashref
  Exceptions : throws if cannot read from $synonyms_file
  Caller     : general
  Status     : Stable

=cut

sub chromosome_synonyms {
  my $self = shift;
  my $file = shift;

  my $synonyms = $self->config->{_chromosome_synonyms} ||= {};

  if ($file) {
    open IN, $file or throw("ERROR: Could not read synonyms file $file: $!");

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
  } elsif (!%{$synonyms} && $self->param('database')) {
    my $sa    = $self->get_adaptor('core', 'slice');
    my $srsa  = $self->get_adaptor('core', 'SeqRegionSynonym');
    $synonyms = _fetch_chr_synonyms($srsa, $sa);
  }
  return $self->config->{_chromosome_synonyms} = $synonyms;
}


=head2 get_source_chr_name
  
  Arg 1      : string $chr
  Arg 2      : (optional) string $set_name
  Arg 3      : (optional) arrayref $valid_chromosomes
  Example    : $syns = $obj->get_source_chr_name()
  Description: Attempts to get the chromosome name as it appears in a source
               (annotation source, FASTA file, etc) from given chromosome name.
               A list of valid chromosome names for the source can be given,
               along with a set name to keep mappings for different sources separate.

               Will also attempt to do simple transforms like adding/removing "chr".

               If no valid match is found, returns $chr as given.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_source_chr_name {
  my ($self, $chr, $set, $valids) = @_;

  $set    ||= 'default';
  $valids ||= [];

  my $chr_name_map = $self->{_chr_name_map}->{$set} ||= {};

  if(!exists($chr_name_map->{$chr})) {
    my $mapped_name = $chr;

    if (!@$valids && $self->can('valid_chromosomes')) {
      my $valid_chr = $self->valid_chromosomes;
      @$valids = ref($valid_chr) eq 'HASH' ? [ keys %$valid_chr ] : @$valid_chr;
    }
    my %valid = map {$_ => 1} @$valids;

    unless($valid{$chr}) {

      # try synonyms first
      my $synonyms = $self->chromosome_synonyms;

      foreach my $syn(keys %{$synonyms->{$chr} || {}}) {
        if($valid{$syn}) {
          $mapped_name = $syn;
          last;
        }
      }

      # still haven't got it
      if($mapped_name eq $chr) {

        # try adding/removing "chr"
        if($chr =~ /^chr/i) {
          my $tmp = $chr;
          $tmp =~ s/^chr//i;

          $mapped_name = $tmp if $valid{$tmp};
        }
        elsif($valid{'chr'.$chr}) {
          $mapped_name = 'chr'.$chr;
        }
      }
    }

    $chr_name_map->{$chr} = $mapped_name;
  }

  return $chr_name_map->{$chr};
}


=head2 add_shortcuts
  
  Arg 1      : string $param_1
  ...
  Arg N      : string $param_N
  Example    : $obj->add_shortcuts($param1, $param2, $param3)
  Description: Adds "shortcuts" to $self for the values of config params
  Returntype : none
  Exceptions : throws if key for $param_N on $self already exists
  Caller     : new() constructors
  Status     : Stable

=cut

sub add_shortcuts {
  my $self = shift;

  return unless $self->config;

  foreach my $param(map {ref($_) eq 'ARRAY' ? @$_ : $_} @_) {
    throw("ERROR: add_shortcuts would overwrite value for \"$param\"\n") if exists($self->{$param});
    $self->{$param} = $self->param($param); 
  }
}


=head2 status_msg
  
  Arg 1      : string $msg
  Example    : $obj->status_msg("Hello world!")
  Description: Prints a time-stamped status message to STDOUT, adding
               a newline character at the end if none given.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub status_msg {
  my $self = shift;
  return if $self->param('quiet');
  
  my $msg = (@_ ? (join "", @_) : "No message");
  print get_time()." - ".$msg.($msg =~ /\n$/ ? "" : "\n");
}


=head2 warning_msg
  
  Arg 1      : string $msg
  Example    : $obj->warning_msg("Look out!")
  Description: Prints a given warning message to handle given
               by $self->warning_fh (default STDERR).
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub warning_msg {
  my $self = shift;
  my $text = shift || '';

  $text = 'WARNING: '.$text unless $text =~ /^warn/i;
  $text = $text."\n" unless $text =~ /\n$/;

  # check if we've seen this warning before
  my $to_check = $text;
  $to_check =~ s/line \d+\n?$//;
  return if $self->{_seen_warnings} && $self->{_seen_warnings}->{$to_check};
  $self->{_seen_warnings}->{$to_check} = 1;

  # $self->config->{warning_count}++ if $self->config;

  my $fh = $self->warning_fh;

  print $fh $text;

  warn($text) unless $self->param('quiet');
}


=head2 skipped_variant_msg

  Arg 1      : string $description
  Arg 2      : (optional) int $line_number
  Example    : $obj->skipped_variant_msg("Cannot parse variant location", 451)
  Description: Prints a given skipped variant as a warning message (based on
               handle given by $self->warning_fh) and writes to file if
               $self->param->{skipped_variants_file} is defined.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub skipped_variant_msg {
  my $self = shift;
  my $description = shift;
  my $line_number = shift || $self->line_number;

  my $msg = $description;
  my $line = $self->{parser}->{current_block} if defined $self->{parser};
  if (defined $line) {
    my $maxChar = 50;
    $line =~ s/\s+/ /g; # trim tabs and extra spaces to improve readability
    $line = substr($line, 0, $maxChar - 4) . "..." if length($line) > $maxChar;
    $msg = "(" . $line . "): " . $msg;
  }

  $msg = (defined $line_number and not defined $self->param('input_data')) ?
    "line $line_number skipped " . $msg : "variant skipped " . $msg;
  $self->warning_msg("WARNING: " . $msg . "\n");

  $self->skipped_variants({
    line => $line,
    line_number => $line_number,
    description => $description
  });
}


=head2 warning_fh
  
  Example    : $fh = $obj->warning_fh()
  Description: Gets filehandle for writing warning messages to.
               Defaults to STDERR, but can be configured via "warning_file" param
  Returntype : glob
  Exceptions : none
  Caller     : warning_msg()
  Status     : Stable

=cut

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


=head2 module_prefix
  
  Example    : $prefix = $obj->module_prefix()
  Description: Gets prefix to be applied when generating classes; this is used
               when instantiating Bio::EnsEMBL::VEP::AnnotationSource::* classes
               so that the correct VEP or Haplo class is created
  Returntype : string
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::AnnotationSourceAdaptor, Bio::EnsEMBL::VEP::CacheDir
  Status     : Stable

=cut

sub module_prefix {
  my $self = shift;
  my $prefix = 'Bio::EnsEMBL::VEP';
  $prefix .= '::Haplo' if $self->param('haplo');
  return $prefix;
}

1;
