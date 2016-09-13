# Copyright [2016] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use Getopt::Long;
use Bio::EnsEMBL::VEP::Haplo::Runner;

my $config = {};

my $arg_count = scalar @ARGV;
my @argv_copy = @ARGV;

GetOptions(
  $config,
  'help',                    # displays help message
  
  # input options,
  'config=s',                # config file name
  'input_file|i=s',          # input file name
  'format=s',                # input file format
  'output_format=s',         # output file format
  'delimiter=s',             # delimiter between fields in input

  # haplo specific
  'haplotype_frequencies=s',
  
  # DB options
  'species=s',               # species e.g. human, homo_sapiens
  'registry=s',              # registry file
  'host=s',                  # database host
  'port=s',                  # database port
  'user|u=s',                  # database user name
  'password=s',              # database password
  'db_version=i',            # Ensembl database version to use e.g. 62
  'assembly|a=s',            # assembly version to use
  'genomes',                 # automatically sets DB params for e!Genomes
  'refseq',                  # use otherfeatures RefSeq DB instead of Ensembl
  'merged',                  # use merged cache
  'all_refseq',              # report consequences on all transcripts in RefSeq cache, includes CCDS, EST etc
  'gencode_basic',           # limit to using just GenCode basic transcript set
  'is_multispecies=i',       # '1' for a multispecies database (e.g protists_euglenozoa1_collection_core_29_82_1)

  # runtime options
  'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
  'phased',                  # force VCF genotypes to be interpreted as phased
  'dont_skip',               # don't skip vars that fail validation
  
  # verbosity options
  'verbose|v',               # print out a bit more info while running
  'quiet',                   # print nothing to STDOUT (unless using -o stdout)
  
  # output options
  'output_file|o=s',         # output file name
  'warning_file=s',          # file to write warnings to
  'force_overwrite',         # force overwrite of output file if already exists
  
  # cache stuff
  'database',                # must specify this to use DB now
  'cache',                   # use cache
  'cache_version=i',         # specify a different cache version
  'show_cache_info',         # print cache info and quit
  'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
  'dir_cache=s',             # specific directory for cache
  'offline',                 # offline mode uses minimal set of modules installed in same dir, no DB connection
  'custom=s' => ($config->{custom} ||= []), # specify custom tabixed bgzipped file with annotationl
  'gff=s',                   # shortcut to --custom [file],,gff
  'gtf=s',                   # shortcut to --custom [file],,gtf
  'fasta=s',                 # file or dir containing FASTA files with reference sequence
  'sereal',                  # user Sereal instead of Storable for the cache
  'synonyms=s',              # file of chromosome synonyms
  
  # debug
  'debug',                   # print out debug info
) or die "ERROR: Failed to parse command-line flags\n";

&usage && exit(0) if (!$arg_count) || $config->{help};

$config->{database} ||= 0;

my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new($config);

$runner->run();



# outputs usage message
sub usage {
    my $usage =<<END;
#-------------#
# HAPLOSAURUS #
#-------------#

version $VERSION
by Will McLaren (wm2\@ebi.ac.uk)

Help: dev\@ensembl.org , helpdesk\@ensembl.org
Twitter: \@ensembl , \@EnsemblWill

Usage:
perl haplo.pl [--cache|--offline|--database] [arguments]

Basic options
=============

--help                 Display this message and quit

-i | --input_file      Input file
-o | --output_file     Output file
--force_overwrite      Force overwriting of output file
--species [species]    Species to use [default: "human"]

END

    print $usage;
}


1;