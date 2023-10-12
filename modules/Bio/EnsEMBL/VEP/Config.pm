=head1 LICENSE

Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

use File::Spec;
use Getopt::Long;
Getopt::Long::Configure("pass_through");

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);


## VALID COMMAND-LINE PARAMETERS IN VEP
#######################################

our @VEP_PARAMS = (
  'help',                    # displays help message

  # input options,
  'config=s',                # config file name
  'input_file|i=s',          # input file name
  'input_data|id=s',         # input data
  'format=s',                # input file format
  'output_format=s',         # output file format
  'delimiter=s',             # delimiter between fields in input
  'no_check_variants_order', # skip check about the variants ordering within a region

  # DB options
  'species|s=s',             # species e.g. human, homo_sapiens
  'registry=s',              # registry file
  'host=s',                  # database host
  'port=s',                  # database port
  'user|u=s',                # database user name
  'password|pass=s',         # database password
  'db_version=i',            # Ensembl database version to use e.g. 62
  'assembly|a=s',            # assembly version to use
  'grch37',                  # set for using GRCh37
  'genomes',                 # automatically sets DB params for e!Genomes
  'refseq',                  # use otherfeatures RefSeq DB instead of Ensembl
  'merged',                  # use merged cache
  'all_refseq',              # report consequences on all transcripts in RefSeq cache, includes CCDS, EST etc
  'gencode_basic',           # limit to using just GenCode basic transcript set
  'is_multispecies=i',       # '1' for a multispecies database (e.g protists_euglenozoa1_collection_core_29_82_1)

  # runtime options
  'transcript_filter=s@',    # filter transcripts
  'exclude_predicted',
  'minimal',                 # convert input alleles to minimal representation
  'most_severe',             # only return most severe consequence
  'summary',                 # only return one line per variation with all consquence types
  'pick',                    # used defined criteria to return most severe line
  'pick_allele',             # choose one con per allele
  'per_gene',                # choose one con per gene
  'pick_allele_gene',        # choose one con per gene, allele
  'flag_pick',               # flag one con per line
  'flag_pick_allele',        # flag one con per allele
  'flag_pick_allele_gene',   # flag one con per gene, allele
  'pick_order=s',            # define the order of categories used by the --*pick* flags
  'buffer_size=i',           # number of variations to read in before analysis
  'failed=i',                # include failed variations when finding existing
  'gp',                      # read coords from GP part of INFO column in VCF (probably only relevant to 1KG)
  'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
  'check_ref',               # check supplied reference allele against DB/FASTA
  'lookup_ref',              # replace supplied reference allele with allele from DB/FASTA
  'check_existing',          # find existing co-located variations
  'check_svs',               # find overlapping structural variations
  'no_check_alleles',        # attribute co-located regardless of alleles
  'exclude_null_alleles',    # exclude variants with null alleles from co-located check (e.g COSMIC)
  'check_frequency',         # enable frequency checking
  'af',                      # add global AF of existing var
  'af_1kg',                  # add 1KG AFs of existing vars
  'af_gnomade',              # add gnomAD v2 exomes AFs of existing vars
  'af_gnomad',               # Same as --af_gnomade to keep old compatibility
  'af_gnomadg',              # add gnomAD v3 genomes AFs of existing vars
  'old_maf',                 # report 1KG/ESP MAFs in the old way (no allele, always < 0.5)
  'max_af',                  # report maximum observed allele frequency in any 1KG, gnomAD v2 exomes, gnomAD v3 genomes pops
  'pubmed',                  # add Pubmed IDs for publications that cite existing vars
  'freq_filter=s',           # exclude or include
  'freq_freq=f',             # frequency to filter on
  'freq_gt_lt=s',            # gt or lt (greater than or less than)
  'freq_pop=s',              # population to filter on
  'filter_common',           # shortcut to MAF filtering
  'allow_non_variant',       # allow non-variant VCF lines through
  'process_ref_homs',        # force processing of individuals with homozygous ref genotype
  'individual=s',            # give results by genotype for individuals
  'individual_zyg=s',        # reports genotypes for individuals
  'phased',                  # force VCF genotypes to be interpreted as phased
  'fork=i',                  # fork into N processes
  'dont_skip',               # don't skip vars that fail validation
  'nearest=s',               # get nearest transcript, gene or symbol (for gene)
  'distance=s',              # set up/downstream distance
  'clin_sig_allele=i',       # use allele specific clinical significance data where it exists
  'overlaps',                # report length and percent of a transcript or regulatory feature overlaped with a SV
  'max_sv_size=i',           # modify the size of structural variant to be handled (limited by default to reduce memory requirements)
  'remove_hgvsp_version',    # removes translation version from hgvs_protein output


  # verbosity options
  'verbose|v',               # print out a bit more info while running
  'quiet|q',                 # print nothing to STDOUT (unless using -o stdout)
  'no_progress',             # don't display progress bars

  # output options
  'everything|e',            # switch on EVERYTHING :-)
  'output_file|o=s@',        # output file name
  'compress_output=s',       # compress output with e.g. bgzip, gzip
  'no_headers',              # don't print headers
  'stats_file|sf=s',         # stats file name
  'stats_text',              # write stats as text
  'stats_html',              # write stats as html
  'no_stats',                # don't write stats file
  'warning_file=s',          # file to write warnings to
  'skipped_variants_file=s', # file name to log skipped variants (not logged otherwise)
  'force_overwrite|force',   # force overwrite of output file if already exists
  'terms|t=s',               # consequence terms to use e.g. NCBI, SO
  'coding_only',             # only return results for consequences in coding regions
  'canonical',               # indicates if transcript is canonical
  'mane',                    # output mane transcript value
  'mane_select',             # output mane select transcript value
  'tsl',                     # output transcript support level
  'appris',                  # output APPRIS transcript annotation
  'ccds',                    # output CCDS identifer
  'xref_refseq',             # output refseq mrna xref
  'uniprot',                 # output Uniprot identifiers (includes UniParc)
  'protein',                 # add e! protein ID to extra column
  'biotype',                 # add biotype of transcript to output
  'hgnc',                    # add HGNC gene ID to extra column
  'symbol',                  # add gene symbol (e.g. HGNC)
  'transcript_version',      # add transcript version to stable id in feature column
  'gene_phenotype',          # indicate if genes are phenotype-associated
  'mirna',                   # identify miRNA structural elements overlapped by variant
  'spdi',                    # add genomic SPDI
  'ga4gh_vrs',               # add GA4GH_VRS
  'hgvs',                    # add HGVS names to extra column
  'hgvsg',                   # add HGVS g. also
  'hgvsg_use_accession',     # force HGVSg to return on chromosome accession instead of input chr name
  'hgvsp_use_prediction',    # force HGVSp to return the notation in predicted format
  'shift_hgvs=i',            # disable/enable 3-prime shifting of HGVS indels to comply with standard
  'ambiguous_hgvs',          # allow input HGVSp. to resolve to many input variants
  'sift=s',                  # SIFT predictions
  'polyphen=s',              # PolyPhen predictions
  'humdiv',                  # use humDiv instead of humVar for PolyPhen
  'condel=s',                # Condel predictions
  'variant_class',           # get SO variant type
  'regulatory',              # enable regulatory stuff
  'cell_type=s',             # filter cell types for regfeats
  'convert=s',               # DEPRECATED: convert input to another format (doesn't run VEP)
  'no_intergenic',           # don't print out INTERGENIC consequences
  'vcf',                     # produce vcf output
  'solr',                    # produce XML output for Solr
  'json',                    # produce JSON document output
  'tab',                     # produce tabulated output
  'vcf_info_field=s',        # allow user to change VCF info field name
  'keep_csq',                # don't nuke existing CSQ fields in VCF
  'keep_ann',                # synonym for keep_csq
  'lrg',                     # enable LRG-based features
  'fields=s',                # define your own output fields
  'domains',                 # output overlapping protein features
  'numbers',                 # include exon and intron numbers
  'total_length',            # give total length alongside positions e.g. 14/203
  'allele_number',           # indicate allele by number to avoid confusion with VCF conversions
  'show_ref_allele',         # indicate reference allele
  'no_escape',               # don't percent-escape HGVS strings
  'ambiguity',               # Add allele ambiguity code
  'var_synonyms', 	         # include variation synonyms in output
  'shift_3prime=i',          # enables shifting of all variants to 3prime
  'shift_genomic=i',         # adds genomic shifting to output, and provides shifting of intergenic variants
  'shift_length',	           # adds the length of the transcript directional shift to output
  'uploaded_allele',         # output allele string given in input 

  # cache stuff
  'database',                # must specify this to use DB now
  'cache:s',                 # use cache (optional param treated like --cache-dir)
  'cache_version=i',         # specify a different cache version
  'show_cache_info',         # print cache info and quit
  'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
  'dir_cache=s',             # specific directory for cache
  'dir_plugins=s',           # specific directory for plugins
  'offline',                 # offline mode uses minimal set of modules installed in same dir, no DB connection
  'fasta|fa=s',              # file or dir containing FASTA files with reference sequence
  'fasta_dir=s',             # dir containing FASTA file (may contain multiple species/assemblies)
  'no_fasta',                # don't autodetect FASTA file in cache dir
  'sereal',                  # user Sereal instead of Storable for the cache
  'synonyms=s',              # file of chromosome synonyms

  # these flags are for use with RefSeq caches
  'bam=s',                   # bam file used to modify transcripts
  'use_transcript_ref',      # extract the reference allele from the transcript (or genome)
  'use_given_ref',           # override use_transcript_ref setting that may be set from cache info

  # custom file stuff
  'custom=s@',               # specify custom tabixed bgzipped or bigWig file with annotation
  'tmpdir=s',                # tmp dir used for BigWig retrieval
  'gff=s',                   # shortcut to --custom [file],,gff
  'gtf=s',                   # shortcut to --custom [file],,gtf
  'bigwig=s',                # shortcut to --custom [file],,bigwig,exact
  'phyloP=s@',               # shortcut to using remote phyloP, may use multiple
  'phastCons=s@',            # shortcut to using remote phastCons, may use multiple
  'ucsc_assembly=s',         # required for phyloP, phastCons, e.g. use hg19 for GRCh37, hg38 for GRCh38
  'ucsc_data_root=s',        # replace if you have the data locally, defaults to http://hgdownload.cse.ucsc.edu/goldenpath/
  'custom_multi_allelic',    # prevents filtering of custom annotation data when comma separated lists are assumed to be allele specific

  # plugins
  'plugin=s@',               # specify a method in a module in the plugins directory
  'safe',                    # die if plugins don't compile or spit warnings

  # debug
  'debug',                   # print out debug info
);


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
  pick_order        => [qw(mane_select mane_plus_clinical canonical appris tsl biotype ccds rank length ensembl refseq )],
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
  individual_zyg
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
      af_gnomade     => 1,
      af_gnomadg     => 1,
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
    flags => ['individual_zyg'],
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
    flags => [qw(check_frequency af af_1kg af_gnomad af_gnomade af_gnomadg max_af pubmed)],
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
      custom => 'file=%gff%,format=gff'
    }
  },
  
  {
    flags => ['gtf'],
    set   => {
      custom => 'file=%gtf%,format=gtf'
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
  pick_order      => [qw(mane_select mane_plus_clinical canonical appris tsl biotype ccds rank length ensembl refseq)],
  nearest         => [qw(transcript gene symbol)],
  compress_output => [qw(gzip bgzip)],
);

# these flags require others to be set, no sensible defaults can be set by us
our %REQUIRES = (
  original  => [qw(filters)],
  phyloP    => [qw(ucsc_assembly)],
  phastCons => [qw(ucsc_assembly)],
  custom_multi_allelic => [qw(custom)]
);

# incompatible options
our %INCOMPATIBLE = (
  most_severe => [qw(biotype no_intergenic protein symbol sift polyphen coding_only ccds mane canonical xref_refseq numbers domains tsl appris uniprot summary pick flag_pick pick_allele flag_pick_allele)],
  summary     => [qw(biotype no_intergenic protein symbol sift polyphen coding_only ccds mane canonical xref_refseq numbers domains tsl appris uniprot most_severe pick flag_pick pick_allele flag_pick_allele)],
  database    => [qw(af_1kg af_gnomad af_gnomade af_gnomadg max_af pubmed var_synonyms offline cache)],
  af_gnomade    => [qw(af_gnomad)],
  quiet       => [qw(verbose)],
  refseq      => [qw(gencode_basic merged)],
  json        => [qw(vcf tab)],
  vcf         => [qw(json tab)],
  tab         => [qw(vcf json)],
  individual  => [qw(minimal individual_zyg)],
  check_ref   => [qw(lookup_ref)],
  check_svs   => [qw(offline)],
  ga4gh_vrs   => [qw(vcf)]
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
  # 4: environment variables
  # 5: %DEFAULTS
  
  # read config file if defined
  $self->read_config_from_file($config->{config}, $config) if defined $config->{config};

  # read default config file if defined
  my $ini_file = ( $config->{dir} || $DEFAULTS{dir} ). '/vep.ini';
  $self->read_config_from_file($ini_file, $config) if -e $ini_file;

  # read environment variables with a VEP_ prefix, e.g. VEP_DIR_PLUGINS
  $self->read_config_from_environment($config);

  # assign default port for GRCh37
  if (defined($config->{'assembly'}) && lc($config->{'assembly'}) eq 'grch37' && defined($config->{'database'}) && !defined($config->{'port'})) {
    $config->{'port'} = 3337;
  }

  # set cache directory based on cache or defaults
  if (defined $config->{cache} && $config->{cache} ne 0){
    $config->{dir_cache} ||= $config->{cache} if -d "$config->{cache}";
    $config->{cache} = 1;
  }

  my $config_command = "";

  my @skip_opts = qw(web_output host port stats_file user warning_file input_data);

  foreach my $flag (sort keys %$config) {
    my $value = $config->{$flag};
    my $default = $DEFAULTS{$flag};
    next if defined($default) && $default eq $value;
    next if !defined($value) || (ref($value) eq "ARRAY" && @{$value} == 0) || grep { /$flag/ } @skip_opts;

    $value = join(" --$flag ", @{$value}) if ref($value) eq "ARRAY";
    
    if ($^O eq "MSWin32"){
      $value =~ s/.+(?=\\)/\[PATH\]/g;
    }
    else {
      $value =~ s/.+(?=\/)/\[PATH\]/g;
    }
    
    $config_command .= $value eq 1? "--$flag "  : "--$flag $value ";
  }

  $config->{full_command} = "vep $config_command";
  $config->{full_command} =~ s/\s*$//;

  # set all other defaults
  foreach my $key(keys %DEFAULTS) {
    $config->{$key} = $DEFAULTS{$key} unless exists($config->{$key});
  }

  # set params that will be arrayrefs empty if not defined
  $config->{$_} ||= [] for @ALLOW_MULTIPLE;
  
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
    delete $config->{$_} for qw(af_1kg af_gnomad af_gnomade af_gnomadg max_af pubmed);
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


=head2 is_valid_param

  Arg 1      : hashref $config
  Arg 2      : string $key
  Arg 3      : string $value
  Example    : $is_valid = is_valid_param($config, "assembly", "GRCh38")
  Description: Check if a param is valid using GetOptions and @VEP_PARAMS, a
               list of VEP's valid commad-line arguments.
  Returntype : bool
  Caller     : read_config_from_environment(), read_config_from_file()
  Status     : Stable

=cut

sub is_valid_param {
  my $config = shift;
  my $key = shift;
  my $value = shift;

  # Prepare variable as command-line arguments
  my @original_ARGV = @ARGV;
  @ARGV = ( "--" . $key, $value );

  # Run GetOptions to check validity of parameter
  my $res = {};
  GetOptions($res, @VEP_PARAMS);
  my $is_valid = %$res ? 1 : 0;

  # Restore command-line arguments
  @ARGV = @original_ARGV;
  return $is_valid;
}


=head2 read_config_from_file

  Arg 1      : string $config_file
  Arg 2      : hashref $config
  Example    : $config_hash = $self->read_config_from_file($config_file, $config)
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
    next if /^\s*$/;

    # preserve spaces between quotes
    s/\s+(?=(?:(?:[^"]*"){2})*[^"]*"[^"]*$)/___SPACE___/g;

    my @split = split /\s+/;
    my $key = shift @split;
    $key =~ s/^\-//g;

    # restore spaces
    s/___SPACE___/ /g for @split;

    # remove quotes
    s/[\"\']//g for @split;

    if(grep {$key eq $_} @ALLOW_MULTIPLE) {
      my $value = join(' ', @split);
      next unless is_valid_param($config, $key, $value);
      push @{$config->{$key}}, $value;
    }
    else {
      next unless is_valid_param($config, $key, $split[0]);
      $config->{$key} ||= $split[0];
    }
  }

  close CONFIG;

  $self->status_msg("Read configuration from $file") if $config->{verbose};

  return $config;
}


=head2 read_config_from_environment

  Arg 1      : hashref $config
  Example    : $config_hash = $self->read_config_from_environment($config)
  Description: Read config params from environment variables prefixed by VEP_
               and add them to config hash; e.g., dir_plugins is set to
               VEP_DIR_PLUGINS, if defined. NB: VEP arguments are assumed to
               always be lowercase.
  Returntype : none
  Caller     : new()
  Status     : Stable

=cut

sub read_config_from_environment {
  my $self = shift;
  my $config = shift;

  for my $var (keys %ENV) {
    # Look for environment variables that start with VEP_
    next unless $var =~ "^VEP_";

    # Avoid setting empty strings
    my $value = $ENV{$var};
    next if $value eq "";

    # Assumption: VEP arguments are always lowercase
    my $key = lc $var;
       $key =~ s/^VEP_//ig;

    my $valid = is_valid_param($config, $key, $value);
    unless ($valid) {
      $self->status_msg("Ignored unsupported option '${key}=${value}' from environment variable $var\n")
        if $config->{verbose};
      next;
    }

    if (grep {$key eq $_} @ALLOW_MULTIPLE) {
      # Properly set flags that can be specified more than once
      push @{$config->{$key}}, $value;

      my $msg = "Appended '${key}=${value}' from environment variable $var (full value of ${key}: ['" .
                join("', '", @{ $config->{$key} }) . "'])\n";
      $self->status_msg($msg) if $config->{verbose};
    } else {
      $config->{$key} ||= $value;
      $self->status_msg("Set '${key}=${value}' from environment variable $var\n")
        if $config->{verbose};
    }
  }

  $self->status_msg("Read configuration from environment variables")
    if $config->{verbose};

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
