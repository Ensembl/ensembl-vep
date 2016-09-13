# Copyright [2016] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Exception;
use FindBin qw($Bin);

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

my $cfg_hash = $test_cfg->base_testing_cfg;
$cfg_hash->{input_file} = $test_cfg->{test_vcf};

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Haplo::Runner');

my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new($cfg_hash);
ok($runner, 'new is defined');

is(ref($runner), 'Bio::EnsEMBL::VEP::Haplo::Runner', 'check class');

is($runner->param('haplo'), 1, 'haplo param set by new');



## METHOD TESTS
###############

is_deeply(
  $runner->get_all_AnnotationSources(),
  [
    bless( {
      '_config' => $runner->config,
      'cache_region_size' => 1000000,
      'dir' => $test_cfg->{cache_dir},
      'serializer_type' => undef,
      'gencode_basic' => undef,
      'source_type' => 'ensembl',
      'all_refseq' => undef,
      'polyphen' => undef,
      'sift' => undef,
      'everything' => undef,
      'info' => {
        'polyphen' => '2.2.2',
        'sift' => 'sift5.2.2',
        'COSMIC' => '75',
        'ESP' => '20141103',
        'gencode' => 'GENCODE 24',
        'HGMD-PUBLIC' => '20154',
        'genebuild' => '2014-07',
        'regbuild' => '13.0',
        'assembly' => 'GRCh38.p5',
        'dbSNP' => '146',
        'ClinVar' => '201601'
      },
      'valid_chromosomes' => [21, 'LRG_485'],
    }, 'Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript' )
  ],
  'get_all_AnnotationSources'
);

# setup_db_connection should return silently in offline mode
ok(!$runner->setup_db_connection(), 'setup_db_connection');

is_deeply($runner->get_valid_chromosomes, [21, 'LRG_485'], 'get_valid_chromosomes');

is_deeply($runner->get_Parser, bless({
  '_config' => $runner->config,
  'file' => $test_cfg->{test_vcf},
  'file_bak' => $test_cfg->{test_vcf},
  'line_number' => 0,
  'check_ref' => undef,
  'chr' => undef,
  'dont_skip' => undef,
  'valid_chromosomes' => {21 => 1, LRG_485 => 1},
  'minimal' => undef,
  'lrg' => undef,
  'delimiter' => "\t",
  'process_ref_homs' => undef,
  'allow_non_variant' => undef,
  'gp' => undef,
  'individual' => undef,
  'phased' => undef,
}, 'Bio::EnsEMBL::VEP::Haplo::Parser::VCF' ), 'get_Parser');

my $t = $runner->get_TranscriptTree;
delete $t->{trees};

is_deeply(
  $runner->get_TranscriptTree,
  bless( {
    '_config' => $runner->config,
    'valid_chromosomes' => {
      '21' => 1,
      'LRG_485' => 1
    },
  }, 'Bio::EnsEMBL::VEP::Haplo::TranscriptTree' ),
  'get_TranscriptTree'
);

is_deeply($runner->get_InputBuffer, bless({
  '_config' => $runner->config,
  'parser' => $runner->get_Parser,
  'buffer_size' => $runner->param('buffer_size'),
  'minimal' => undef,
  'transcript_tree' => $runner->get_TranscriptTree,
}, 'Bio::EnsEMBL::Haplo::VEP::InputBuffer' ), 'get_InputBuffer');


throws_ok {$runner->haplotype_frequencies($test_cfg->{test_gzvcf})} qr/No hex header found/, 'haplotype_frequencies - no hex header found';
throws_ok {$runner->haplotype_frequencies($test_cfg->{chr_synonyms})} qr/first column should be an md5sum/, 'haplotype_frequencies - no hex column found';

ok($runner->haplotype_frequencies($test_cfg->{haplo_freqs}), 'haplotype_frequencies');

is_deeply(
  $runner->haplotype_frequencies->{'79d615f301917d96d4cefaea29702508'},
  [
    {
      transcript => 'ENST00000307301',
      ALL => 0.889,
      AFR => 0.708,
      AMR => 0.952,
      EAS => 0.944,
      EUR => 0.946,
      SAS => 0.971
    }
  ],
  'haplotype_frequencies - check data'
);


ok($runner->init, 'init');


## run method
$runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({%$cfg_hash, output_file => 'STDOUT'});

no warnings 'once';
open(SAVE, ">&STDOUT") or die "Can't save STDOUT\n"; 

my $tmp;
close STDOUT;
open STDOUT, '>', \$tmp;

$runner->run();

$tmp =~ s/\t/  /g;
is(
  $tmp,
qq{ENST00000352957  ENST00000352957:612A>G    ENSP00000284967:REF      HG00096  1
ENST00000352957  ENST00000352957:91T>C,612A>G    ENSP00000284967:31S>P      HG00096  1
ENST00000307301  ENST00000307301:612A>G    ENSP00000305682:REF      HG00096  1
ENST00000307301  ENST00000307301:91T>C,612A>G    ENSP00000305682:31S>P      HG00096  1
ENST00000419219  ENST00000419219:582A>G    ENSP00000404426:REF      HG00096  1
ENST00000419219  ENST00000419219:91T>C,582A>G    ENSP00000404426:31S>P      HG00096  1
},
  'run - check output'
);


$runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({%$cfg_hash, output_file => 'STDOUT'});
$tmp = undef;

$runner->haplotype_frequencies($test_cfg->{haplo_freqs});
$runner->run;

is_deeply(
  {
    map {(split("=", $_))[0] => (split("=", $_))[1]}
    map {split(",", $_)}
    map {(split("\t", $_))[5]}
    grep {/ENST00000307301\:612/} 
    split("\n", $tmp)
  },
  {
    'ALL' => '0.0905',
    'AMR' => '0.0432',
    'AFR' => '0.234',
    'EUR' => '0.0467',
    'SAS' => '0.0194',
    'EAS' => '0.0476'
  },
  'run with haplotype_frequencies'
);


# restore STDOUT
open(STDOUT, ">&SAVE") or die "Can't restore STDOUT\n";


done_testing();
