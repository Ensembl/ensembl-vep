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

# EnsEMBL module for Bio::EnsEMBL::VEP::Config
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Config - Class used to configure VEP

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Config;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

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
  pick_order        => [qw(canonical appris tsl biotype ccds rank length ensembl refseq)],
  terminal_width    => 48,
  vcf_info_field    => 'CSQ',
  
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
);

# these flags can be specified more than once on the command line
# and so get turned into listrefs
our @ALLOW_MULTIPLE = qw(
  plugin
  custom
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
      af_exac        => 1,
      pubmed         => 1,
      uniprot        => 1,
      tsl            => 1,
      appris         => 1,
      variant_class  => 1,
      gene_phenotype => 1,
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
    flags => ['check_frequency'],
    set => {
      no_check_alleles => 1,
    }
  },
  
  {
    flags => [qw(check_frequency af af_1kg af_esp af_exac pubmed)],
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
    flags => ['gvf'],
    set   => {
      output_format => 'gvf',
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
);

# valid values for certain flags
our %VALID = (
  format     => [qw(ensembl vcf hgvs id pileup guess)],
  convert    => [qw(ensembl vcf hgvs pileup)],
  terms      => [qw(SO display NCBI)],
  sift       => [qw(s p b)],
  polyphen   => [qw(s p b)],
  pick_order => [qw(canonical appris tsl biotype ccds rank length ensembl refseq)],
);

# incompatible options
our %INCOMPATIBLE = (
  most_severe => [qw(no_intergenic protein symbol sift polyphen coding_only ccds canonical xref_refseq numbers domains tsl appris uniprot summary pick flag_pick pick_allele flag_pick_allele)],
  summary     => [qw(no_intergenic protein symbol sift polyphen coding_only ccds canonical xref_refseq numbers domains tsl appris uniprot most_severe pick flag_pick pick_allele flag_pick_allele)],
  database    => [qw(af_1kg af_esp af_exac pubmed offline cache)],
  quiet       => [qw(verbose)],
  refseq      => [qw(gencode_basic merged)],
  merged      => [qw(database)],
  json        => [qw(vcf gvf tab)],
  vcf         => [qw(json gvf tab)],
  gvf         => [qw(vcf json tab)],
  tab         => [qw(vcf gvf json)],
);

# deprecated/replaced flags
our %DEPRECATED = (
  'gmaf'     => 'af',
  'maf_1kg'  => 'af_1kg',
  'maf_esp'  => 'af_1kg',
  'maf_exac' => 'af_exac',
  'check_alleles' => undef,
  'html'     => undef,
  'gvf'      => undef,
  'convert'  => undef,
);


####################################
####################################
####################################



## METHODS
##########

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
  foreach my $key(keys %DEFAULTS) {
    $config->{$key} = $DEFAULTS{$key} unless exists($config->{$key});
  }
    
  # config file?
  if(defined $config->{config}) {
    my $from_file = $self->read_config_from_file($config->{config});

    # copy to working config hashref
    $config->{$_} = $from_file->{$_} for keys %$from_file;
  }
  
  # ini file?
  my $ini_file = $config->{dir}.'/vep.ini';

  if(-e $ini_file) {
    my $from_file = $self->read_config_from_file($ini_file);

    # copy to working config hashref
    $config->{$_} = $from_file->{$_} for keys %$from_file;
  }
  
  # these flags need turning into listrefs
  foreach my $flag(grep {defined($config->{$_}) && ref($config->{$_}) ne 'ARRAY'} @LIST_FLAGS) {
    $config->{$flag} = [split(',', $config->{$flag})];
  }
  
  # apply option sets
  foreach my $set(@OPTION_SETS) {
    foreach my $flag(@{$set->{flags}}) {
      if(
        defined($config->{$flag}) &&
        (
          ref($config->{$flag}) eq 'ARRAY' ?
          scalar @{$config->{$flag}} :
          $config->{$flag}
        )
      ) {
        foreach my $key(keys %{$set->{set}}) {
          my $val = $set->{set}->{$key};
          $val =~ s/\%(\S+)\%/$config->{$1}/g;
          
          if(ref($config->{$key}) eq 'ARRAY') {
            push @{$config->{$key}}, $val;
          }
          else {
            $config->{$key} = $val;
          }
        }
      }
    }
  }
  
  # check config before we return
  $self->check_config($config);

  # copy to $self
  $self->{_params} = $config;
    
  return $self;
}

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
    delete $config->{$_} for qw(af_1kg af_esp af_exac pubmed);
  }
  
  # check valid values for flags
  foreach my $flag(grep {defined($config->{$_})} keys %VALID) {
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
  unless($config->{safe}) {
    foreach my $flag(grep {$config->{$_}} keys %INCOMPATIBLE) {
      foreach my $invalid(grep {$config->{$_}} @{$INCOMPATIBLE{$flag}}) {
        throw sprintf("ERROR: Can't use --%s and --%s together\n", $flag, $invalid);
      }
    }
  }

  # check deprecated flags
  foreach my $flag(grep {$config->{$_}} keys %DEPRECATED) {
    my $msg = "ERROR: --$flag has been deprecated";
    
    if(my $new = $DEPRECATED{$flag}) {
      $msg .= " - please use --$new instead";
    }

    throw("$msg\n");
  }

  ## HACK FIX FOR UNEXPLAINED WEB BUG
  # sometimes web jobs e.g. RT 171600 see "Cache directory [blah]/homo_sapiens1_merged not found" errors
  # "1" is getting appended to the species name somehow
  $config->{species} =~ s/\d+$//;
  
  # check one of database/cache/offline/custom
  if(!grep {$config->{$_}} qw(database cache offline custom)) {
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
Installer: http://www.ensembl.org/info/docs/tools/vep/vep_script.html#installer
Cache: http://www.ensembl.org/info/docs/tools/vep/script/index.html#cache

    }
  };
  
  # check if using filter and original
  throw("ERROR: You must also provide output filters using --filter to use --original\n") if defined($config->{original}) && !defined($config->{filter});
}

# reads config from a file
sub read_config_from_file {
  my $self = shift;
  my $file = shift;

  throw("ERROR: Could not open config file \"".($file || '')."\"\n") unless $file;

  open CONFIG, $file or throw("ERROR: Could not open config file \"$file\"\n");

  my $config = {};

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

# gets/sets the value of a config parameter given a key
sub param {
  my $self = shift;
  my $key = shift;

  throw("No parameter name given") unless $key;

  $self->{_params}->{$key} = shift if @_;

  return $self->{_params}->{$key};
}

# returns raw unaltered config hash as given by user
# this will typically be the hash derived from GetOpt
# in the VEP script
sub _raw_config {
  return $_[0]->{_raw_config};
}

1;
