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

# EnsEMBL module for Bio::EnsEMBL::VEP::Config
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Config - Class used to configure VEP

=head1 SYNOPSIS

my $config = Bio::EnsEMBL::VEP::Config->new();

my $species = $config->param('species');

=head1 DESCRIPTION

The Config class is used to store default and user parameters required
for configuring a VEP runner. A reference to a single instance is passed
around all classes derived from a VEP runner.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Config;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);


## GLOBAL VARIABLES USED FOR INITIAL CONFIG AND SETUP
#####################################################

# default settings
our %DEFAULTS = (  

  # database settings
  database          => 1,
  host              => 'ensembldb.ensembl.org',
  user              => 'anonymous',
  port              => 3306,
  species           => 'homo_sapiens',
  no_slice_cache    => 0,
  
  # cache settings
  dir               => $ENV{HOME}.'/.vep',
  chunk_size        => 50000,
  cache_region_size => 1000000,
  
  # misc settings
  distance          => $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE.','.$Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE,
  buffer_size       => 5000,
  input_file        => 'STDIN',
  delimiter         => "\t",
  output_file       => "variant_effect_output.txt",
  tmpdir            => '/tmp',
  format            => 'guess',
  output_format     => 'vep',
  terms             => 'SO',
  failed            => 0,
  core_type         => 'core',
  polyphen_analysis => 'humvar',
  pick_order        => [qw(canonical appris tsl biotype ccds rank length ensembl refseq mane)],
  terminal_width    => 48,
  vcf_info_field    => 'CSQ',
  ucsc_data_root    => 'http://hgdownload.cse.ucsc.edu/goldenpath/',
  max_sv_size       => 10000000,
  clin_sig_allele   => 1,
  
  # frequency filtering
  freq_freq         => 0.01,
  freq_filter       => 'exclude',
  freq_pop          => '1KG_ALL',
  freq_gt_lt        => 'gt',
);

# these flags take comma-separated lists
# and so need to be converted to listrefs
our @LIST_FLAGS = qw(
  individual
  cell_type
  pick_order
  fields
  chr
  distance
  dont_export
);

# these flags can be specified more than once on the command line
# and so get turned into listrefs
our @ALLOW_MULTIPLE = qw(
  plugin
  custom
  phyloP
  phastCons
  transcript_filter
);

# sets of options that turn on / off others
# if user sets any of the flags in "flags", all the flags in "set" get applied
# NB: order is important, hence why this is a list!
our @OPTION_SETS = (
  {
    flags => ['everything'],
    set   => {
      sift           => 'b',
      polyphen       => 'b',
      ccds           => 1,
      hgvs           => 1,
      symbol         => 1,
      numbers        => 1,
      domains        => 1,
      regulatory     => 1,
      canonical      => 1,
      protein        => 1,
      biotype        => 1,
      af             => 1,
      af_1kg         => 1,
      af_esp         => 1,
      af_gnomad      => 1,
      max_af         => 1,
      pubmed         => 1,
      uniprot        => 1,
      mane           => 1,
      tsl            => 1,
      appris         => 1,
      variant_class  => 1,
      gene_phenotype => 1,
      mirna          => 1,
    },
  },
  
  {
    flags => ['hgvs'],
    set   => {
      hgvsc => 1,
      hgvsp => 1,
    },
  },
  
  {
    flags => ['genomes'],
    set   => {
      host => 'mysql-eg-publicsql.ebi.ac.uk',
      port => 4157,
    },
  },
  
  {
    flags => ['grch37'],
    set   => {
      port => 3337,
      assembly => 'GRCh37',
    },
  },
  
  {
    flags => ['refseq'],
    set   => {
      core_type => 'otherfeatures',
    },
  },
  
  {
    flags => ['individual'],
    set   => {
      allow_non_variant => 1,
    },
  },
  
  {
    flags => ['cell_type'],
    set   => {
      regulatory => 1,
    },
  },
  
  {
    flags => ['filter_common'],
    set   => {
      check_frequency => 1,
    },
  },
  
  {
    flags => [qw(check_frequency af af_1kg af_esp af_exac af_gnomad max_af pubmed)],
    set   => {
      check_existing => 1,
    },
  },

  {
    flags => [qw(show_cache_info)],
    set   => {
      offline => 1,
      format => 'ensembl',
    }
  },
  
  {
    flags => [qw(offline)],
    set   => {
      cache => 1,
      database => 0,
    }
  },
  
  {
    flags => ['cache'],
    set   => {
      no_slice_cache => 1,
    }
  },
  
  {
    flags => ['convert'],
    set   => {
      no_stats => 1,
    }
  },
  
  {
    flags => ['json'],
    set   => {
      rest => 1,
      output_format => 'json',
    }
  },
  
  {
    flags => ['rest'],
    set   => {
      no_escape => 1,
    }
  },
  
  {
    flags => ['vcf'],
    set   => {
      output_format => 'vcf',
      symbol => 1,
      biotype => 1,
      numbers => 1
    }
  },
  
  {
    flags => ['tab'],
    set   => {
      output_format => 'tab',
    }
  },
  
  {
    flags => ['humdiv'],
    set   => {
      polyphen_analysis => 'humdiv',
    }
  },
  
  {
    flags => ['minimal'],
    set   => {
      allele_number => 1,
    }
  },
  
  {
    flags => ['gff'],
    set   => {
      custom => '%gff%,,gff'
    }
  },
  
  {
    flags => ['gtf'],
    set   => {
      custom => '%gtf%,,gtf'
    }
  },
  
  {
    flags => ['bigwig'],
    set   => {
      custom => '%bigwig%,,bigwig,exact'
    }
  },
  
  {
    flags => ['phyloP'],
    set   => {
      custom => '%ucsc_data_root%%ucsc_assembly%/phyloP%phyloP%way/%ucsc_assembly%.phyloP%phyloP%way.bw,phlyoP%phyloP%way,bigwig,exact'
    }
  },
  
  {
    flags => ['phastCons'],
    set   => {
      custom => '%ucsc_data_root%%ucsc_assembly%/phastCons%phastCons%way/%ucsc_assembly%.phastCons%phastCons%way.bw,phastCons%phastCons%way,bigwig,exact'
    }
  },

  {
    flags => ['exclude_predicted'],
    set   => {
      transcript_filter => 'not stable_id match ^X._'
    }
  },

  {
    flags => ['bam'],
    set   => {
      use_transcript_ref => 1,
      bam_edited => 1,
    }
  },
  {
    flags => ['mane'],
    set   => {
      mane_select => 1,
      mane_plus_clinical => 1,
    }
  },
  {
    flags => ['use_given_ref'],
    set   => {
      use_transcript_ref => 0,
    }
  },
);

# valid values for certain flags
our %VALID = (
  format          => [qw(ensembl vcf hgvs id spdi region guess)],
  terms           => [qw(SO display NCBI)],
  sift            => [qw(s p b)],
  polyphen        => [qw(s p b)],
  pick_order      => [qw(mane canonical appris tsl biotype ccds rank length ensembl refseq)],
  nearest         => [qw(transcript gene symbol)],
  compress_output => [qw(gzip bgzip)],
);

# these flags require others to be set, no sensible defaults can be set by us
our %REQUIRES = (
  original  => [qw(filters)],
  phyloP    => [qw(ucsc_assembly)],
  phastCons => [qw(ucsc_assembly)],
  custom_multi_allelic => [qw(custom)],
);

# incompatible options
our %INCOMPATIBLE = (
  most_severe => [qw(biotype no_intergenic protein symbol sift polyphen coding_only ccds mane canonical xref_refseq numbers domains tsl appris uniprot summary pick flag_pick pick_allele flag_pick_allele)],
  summary     => [qw(biotype no_intergenic protein symbol sift polyphen coding_only ccds mane canonical xref_refseq numbers domains tsl appris uniprot most_severe pick flag_pick pick_allele flag_pick_allele)],
  database    => [qw(af_1kg af_esp af_exac af_gnomad max_af pubmed var_synonyms offline cache)],
  quiet       => [qw(verbose)],
  refseq      => [qw(gencode_basic merged)],
  json        => [qw(vcf tab)],
  vcf         => [qw(json tab)],
  tab         => [qw(vcf json)],
  individual  => [qw(minimal)],
  check_ref   => [qw(lookup_ref)],
  check_svs   => [qw(offline)],
);

# deprecated/replaced flags
# e.g. 'gmaf' => 'af'
our %DEPRECATED = (
  'convert' => undef,
);

our %UNSORTABLE = (
  id     => 1,
  hgvs   => 1,
  spdi   => 1,
  region => 1,
);

####################################
####################################
####################################



## METHODS
##########


=head2 new

  Arg 1      : hashref $args
               {
                 arg1 => val1,
                 ...
                 argN => valN,
               }
  Example    : $config = Bio::EnsEMBL::VEP::Config->new($args);
  Description: Constructor for Bio::EnsEMBL::VEP::Config. The args hashref
               may contain any parameters that you wish to add or modify.
               Sensible defaults should be set for most if not all parameters,
               see the %DEFAULTS hash.
  Returntype : Bio::EnsEMBL::VEP::Config
  Exceptions : throws if first arg is not hashref
  Caller     : Runner constructor
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # initialise self
  # we will usually be passed a hashref of options
  # this will typically be from GetOpt
  my $self = bless {_raw_config => shift || {}}, $class;

  throw('First argument given is not a hashref') unless ref($self->_raw_config) eq 'HASH';
  
  # make a hash copy of the raw
  # we might need the raw initial config later?
  my $config;
  %$config = %{$self->_raw_config};
  
  ## FLAG PRIORITY
  # 1: command line flags
  # 2: config file
  # 3: $HOME/.vep/vep.ini file
  # 4: %DEFAULTS
  
  # set defaults
  # we need to do this first so paths are set for reading config etc

  # Assign default port for GRCh37
  if (defined($config->{'assembly'}) && lc($config->{'assembly'}) eq 'grch37' && defined($config->{'database'}) && !defined($config->{'port'})) {
    $config->{'port'} = 3337;
  }

  foreach my $key(keys %DEFAULTS) {
    $config->{$key} = $DEFAULTS{$key} unless exists($config->{$key});
  }

  # set those that will be arrayrefs empty if not defined
  $config->{$_} ||= [] for @ALLOW_MULTIPLE;
    
  # config file?
  $self->read_config_from_file($config->{config}, $config) if defined $config->{config};
  
  # ini file?
  my $ini_file = $config->{dir}.'/vep.ini';
  $self->read_config_from_file($ini_file, $config) if -e $ini_file;
  
  # these flags need turning into listrefs
  foreach my $flag(grep {defined($config->{$_}) && ref($config->{$_}) ne 'ARRAY'} @LIST_FLAGS) {
    $config->{$flag} = [split(',', $config->{$flag})];
  }

  $config->{unsorted_formats} = \%UNSORTABLE;

  $self->apply_option_sets($config);
  
  # check config before we return
  $self->check_config($config);

  # copy to $self
  $self->{_params} = $config;
    
  return $self;
}


=head2 apply_option_sets

  Arg 1      : hashref $config_hashref
  Example    : $config->apply_option_sets($config_hashref)
  Description: Applies "option sets". These allow code to switch
               on or set the value of one or more flags if a primary
               flag is set. Option sets are hard-coded in @OPTION_SETS.
  Returntype : hashref $config_hashref
  Exceptions : none
  Caller     : new()
  Status     : Stable

=cut

sub apply_option_sets {
  my $self = shift;
  my $config = shift;

  # apply option sets
  foreach my $set(@OPTION_SETS) {
    foreach my $flag(@{$set->{flags}}) {

      # check whether the flag has been set
      if($self->_is_flag_active($config, $flag)) {

        # some flags are allowed to be set more than once
        # so we're going to iterate over each instance specified
        my @instances = (
          ref($config->{$flag}) eq 'ARRAY' ?
          @{$config->{$flag}} :
          ($config->{$flag})
        );

        foreach my $i(0..$#instances) {

          # now we work on the set of flags that are going to be altered
          foreach my $key(keys %{$set->{set}}) {
            my $val = $set->{set}->{$key};

            # the syntax allows us to replace values with other config params
            # for example, 'blah_%foo%_blah' will have '%foo%' replaced by $config->{foo}
            my %replace;

            # again though, we have to take care of the cases where the same flag can be set multiple times
            while($val =~ m/\%(\S+?)\%/g) {
              my $sub;
              next unless defined($config->{$1});

              # for those set multiple times
              if(ref($config->{$1}) eq 'ARRAY') {

                # replace from the appropriate instance of our "foo"
                # defaulting to the first if further ones are not defined
                $sub = defined($config->{$1}->[$i]) ? $config->{$1}->[$i] : $config->{$1}->[0];
              }

              # for those set once
              else {
                $sub = $config->{$1};
              }

              $replace{$1} = $sub;
            }

            # now do the actual replacement
            $val =~ s/\%(\S+?)\%/$replace{$1}||"%$1%"/ge;
            
            # and set the value, either by appending to the arrayref
            if(ref($config->{$key}) eq 'ARRAY') {
              push @{$config->{$key}}, $val;
            }
            # or setting the value in $config
            else {
              $config->{$key} = $val;
            }
          }
        }
      }
    }
  }

  return $config
}


=head2 _is_flag_active

  Arg 1      : hashref $config_hashref
  Arg 2      : string $flag
  Example    : $config->apply_option_sets($config_hashref)
  Description: Resolve whether a flag has been set. This depends on
               if the value is defined and if defined whether it is
               an arrayref (empty arrayref counts as unset).
  Returntype : bool
  Exceptions : none
  Caller     : check_config(), apply_option_sets()
  Status     : Stable

=cut

sub _is_flag_active {
  my ($self, $config, $flag) = @_;
  return 0 unless defined($config->{$flag});
  return ref($config->{$flag}) eq 'ARRAY' ? scalar @{$config->{$flag}} : $config->{$flag};
}


=head2 check_config

  Arg 1      : hashref $config_hashref
  Example    : $config->check_config($config_hashref)
  Description: Applies various checks to a config hashref before
               constructor returns successfully:
                - turn on quiet mode when using STDOUT as output
                - disable certain options if using --everything and --database
                - check for valid values for flags using %VALID
                - check for incompatible flags using %INCOMPATIBLE
                - check for required flags using %REQUIRED
                - check for deprecated flags using %DEPRECATED
                - check if one of --database, --cache, --offline or --custom is specified
  Returntype : none
  Exceptions : throws if any of above checks fail
  Caller     : new()
  Status     : Stable

=cut

sub check_config {
  my $self = shift;
  my $config = shift;

  # force quiet if outputting to STDOUT
  if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
    delete $config->{verbose} if defined($config->{verbose});
    $config->{quiet} = 1;
  }
  
  # turn off some options if using --everything and --database
  if($config->{everything} && $config->{database}) {
    delete $config->{$_} for qw(af_1kg af_esp af_exac af_gnomad max_af pubmed);
  }
  
  # check valid values for flags
  foreach my $flag(grep {$self->_is_flag_active($config, $_)} keys %VALID) {
    my @values = ref($config->{$flag}) eq 'ARRAY' ? @{$config->{$flag}} : split(',', $config->{$flag});
    
    foreach my $value(@values) {
      throw
        sprintf(
          "ERROR: \"%s\" is not a valid value for --%s\nValid values: %s\n",
          $value,
          $flag,
          join(", ", @{$VALID{$flag}})
        ) unless grep {$value eq $_} @{$VALID{$flag}};
    }
  }
  
  # check incompatible flags
  # exception: var_synonyms works online only for Variant Recoder
  unless($config->{safe}) {
    foreach my $flag(grep {$self->_is_flag_active($config, $_)} keys %INCOMPATIBLE) {
      foreach my $invalid(grep {$self->_is_flag_active($config, $_)} @{$INCOMPATIBLE{$flag}}) {
        die sprintf("ERROR: Can't use --%s and --%s together\n", $flag, $invalid) unless $self->{_raw_config}->{is_vr} && $flag eq "database" && $invalid eq "var_synonyms";
      }
    }
  }
  
  # check required flags
  foreach my $flag(grep {$self->_is_flag_active($config, $_)} keys %REQUIRES) {
    foreach my $required(@{$REQUIRES{$flag}}) {
      die sprintf("ERROR: You must set --%s to use --%s\n", $required, $flag) unless $self->_is_flag_active($config, $required);
    }
  }

  # check deprecated flags
  foreach my $flag(grep {$self->_is_flag_active($config, $_)} keys %DEPRECATED) {
    my $msg = "ERROR: --$flag has been deprecated";
    
    if(my $new = $DEPRECATED{$flag}) {
      $msg .= " - please use --$new instead";
    }

    die "$msg\n";
  }

  ## HACK FIX FOR UNEXPLAINED WEB BUG
  # sometimes web jobs e.g. RT 171600 see "Cache directory [blah]/homo_sapiens1_merged not found" errors
  # "1" is getting appended to the species name somehow
  # $config->{species} =~ s/\d+$//;
  
  # check one of database/cache/offline/custom
  if(!grep {$self->_is_flag_active($config, $_)} qw(database cache offline custom)) {
    die qq{
IMPORTANT INFORMATION:

The VEP can read gene data from either a local cache or local/remote databases.

Using a cache is the fastest and most efficient way to use the VEP. The
included INSTALL.pl script can be used to fetch and set up cache files from the
Ensembl FTP server. Simply run "perl INSTALL.pl" and follow the instructions, or
see the documentation pages listed below.

If you have already set up a cache, use "--cache" or "--offline" to use it.

It is possible to use the public databases hosted at ensembldb.ensembl.org, but
this is slower than using the cache and concurrent and/or long running VEP jobs
can put strain on the Ensembl servers, limiting availability to other users.

To enable using databases, add the flag "--database".

Documentation
Installer: http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer
Cache: http://www.ensembl.org/info/docs/tools/vep/script/index.html#cache

    }
  };
}


=head2 read_config_from_file

  Arg 1      : string $config_file
  Arg 2      : hashref $config
  Example    : $config_hash = $config->read_config_from_file($config_file, $config)
  Description: Read config params from a flat file and add them to config hash
  Returntype : none
  Exceptions : throws if cannot read from file
  Caller     : new()
  Status     : Stable

=cut

sub read_config_from_file {
  my $self = shift;
  my $file = shift;
  my $config = shift;

  throw("ERROR: Could not open config file \"".($file || '')."\"\n") unless $file;

  open CONFIG, $file or throw("ERROR: Could not open config file \"$file\"\n");

  while(<CONFIG>) {
    next if /^\#/;

    # preserve spaces between quotes
    s/([\"\'].*)(\s)(.*[\"\'])/$1\_\_\_SPACE\_\_\_$3/g;

    my @split = split /\s+|\=/;
    my $key = shift @split;
    $key =~ s/^\-//g;

    # restore spaces
    s/\_\_\_SPACE\_\_\_/ /g for @split;

    # remove quotes
    s/[\"\']//g for @split;

    if(grep {$key eq $_} @ALLOW_MULTIPLE) {
      push @{$config->{$key}}, @split;
    }
    else {
      $config->{$key} ||= $split[0];
    }
  }

  close CONFIG;

  $self->status_msg("Read configuration from $file") if $self->param('verbose');

  return $config;
}


=head2 param

  Arg 1      : string $param_name
  Arg 2      : (optional) $new_value
  Example    : $value = $config->param($param)
               $config->param($param, $new_value)
  Description: Getter/setter for config parameters
  Returntype : none
  Exceptions : throws if no parameter name given
  Caller     : BaseVEP param()
  Status     : Stable

=cut

sub param {
  my $self = shift;
  my $key = shift;

  throw("No parameter name given") unless $key;

  $self->{_params}->{$key} = shift if @_;

  return $self->{_params}->{$key};
}


=head2 _raw_config

  Example    : $raw = $config->_raw_config();
  Description: Gets the raw config hash as originally supplied by the user
               to the new() method
  Returntype : hashref
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _raw_config {
  return $_[0]->{_raw_config};
}

1;
