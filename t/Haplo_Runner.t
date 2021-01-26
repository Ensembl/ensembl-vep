# Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
use Test::Warnings qw(warning :no_end_test);
use FindBin qw($Bin);

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

my $cfg_hash = $test_cfg->base_testing_cfg;
$cfg_hash->{input_file} = $test_cfg->{test_vcf};

## BASIC TESTS
##############

use_ok('Bio::EnsEMBL::VEP::AnnotationType::Transcript');

SKIP: {

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  no warnings 'once';
  skip 'Set::IntervalTree not installed', 26 unless $Bio::EnsEMBL::VEP::AnnotationType::Transcript::CAN_USE_INTERVAL_TREE;

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
        'filter' => [],
        'info' => {
          'sift' => 'sift5.2.2',
          'polyphen' => '2.2.2',
          '1000genomes' => 'phase3',
          'COSMIC' => '80',
          'ESP' => 'V2-SSA137',
          'gnomAD' => '170228',
          'gencode' => 'GENCODE 24',
          'genebuild' => '2014-07',
          'HGMD-PUBLIC' => '20164',
          'regbuild' => '13.0',
          'dbSNP' => '149',
          'ClinVar' => '201704',
          'assembly' => 'GRCh38.p5'
        },
        'valid_chromosomes' => [21, 22, 'LRG_485'],
        'bam' => undef,
        'use_transcript_ref' => undef,
        'nearest' => undef,
      }, 'Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript' )
    ],
    'get_all_AnnotationSources'
  );

  # setup_db_connection should return silently in offline mode
  ok(!$runner->setup_db_connection(), 'setup_db_connection');

  is_deeply($runner->valid_chromosomes, [21, 22, 'LRG_485'], 'valid_chromosomes');

  is_deeply($runner->get_Parser, bless({
    '_config' => $runner->config,
    'file' => $test_cfg->{test_vcf},
    'file_bak' => $test_cfg->{test_vcf},
    'line_number' => 0,
    'check_ref' => undef,
    'lookup_ref' => undef,
    'chr' => undef,
    'dont_skip' => undef,
    'valid_chromosomes' => {21 => 1, 22 => 1, LRG_485 => 1},
    'minimal' => undef,
    'lrg' => undef,
    'delimiter' => "\t",
    'process_ref_homs' => undef,
    'allow_non_variant' => undef,
    'gp' => undef,
    'individual' => undef,
    'phased' => undef,
     'max_sv_size' => 10000000,
  }, 'Bio::EnsEMBL::VEP::Haplo::Parser::VCF' ), 'get_Parser');

  my $t = $runner->get_TranscriptTree;
  delete $t->{trees};

  is_deeply(
    $runner->get_TranscriptTree,
    bless( {
      '_config' => $runner->config,
      'valid_chromosomes' => [21, 22, 'LRG_485'],
    }, 'Bio::EnsEMBL::VEP::TranscriptTree' ),
    'get_TranscriptTree'
  );

  is_deeply($runner->get_InputBuffer, bless({
    '_config' => $runner->config,
    'parser' => $runner->get_Parser,
    'buffer_size' => $runner->param('buffer_size'),
    'minimal' => undef,
    'max_not_ordered_variants' => 100,
    'max_not_ordered_variants_distance' => 100,
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
    'haplotype_frequencies - check data 1'
  );

  is_deeply(
    $runner->haplotype_frequencies->{'91212e92e2a1940184a2eb77fd115d5c'},
    [
      {
        transcript => 'ENST00000307301',
        ALL => '0.0905',
        AMR => '0.0432',
        AFR => '0.234',
        EUR => '0.0467',
        SAS => '0.0194',
        EAS => '0.0476'
      },
      {
        transcript => 'ENSTtest',
        ALL => '0.0905',
        AMR => '0.0432',
        AFR => '0.234',
        EUR => '0.0467',
        SAS => '0.0194',
        EAS => '0.0476'
      }
    ],
    'haplotype_frequencies - check data 2'
  );

  $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new($cfg_hash);
  ok($runner->haplotype_frequencies($test_cfg->{haplo_freqs_nh}), 'haplotype_frequencies - no header');

  is_deeply(
    $runner->haplotype_frequencies->{'79d615f301917d96d4cefaea29702508'},
    [
      {
        unknown1 => 'ENST00000307301',
        freq1 => 0.889,
        freq2 => 0.708,
        freq3 => 0.952,
        freq4 => 0.944,
        freq5 => 0.946,
        freq6 => 0.971
      }
    ],
    'haplotype_frequencies - no header - check data'
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
  qq{ENST00000352957  ENST00000352957:612A>G    ENSP00000284967:REF      rs1135618  HG00096:1
ENST00000352957  ENST00000352957:91T>C,612A>G    ENSP00000284967:31S>P      rs3989369,rs1135618  HG00096:1
ENST00000307301  ENST00000307301:612A>G    ENSP00000305682:REF      rs1135618  HG00096:1
ENST00000307301  ENST00000307301:91T>C,612A>G    ENSP00000305682:31S>P      rs3989369,rs1135618  HG00096:1
ENST00000419219  ENST00000419219:582A>G    ENSP00000404426:REF      rs1135618  HG00096:1
ENST00000419219  ENST00000419219:91T>C,582A>G    ENSP00000404426:31S>P      rs3989369,rs1135618  HG00096:1
},
    'run - check output'
  );


  $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({%$cfg_hash, output_file => 'STDOUT'});
  $tmp = undef;
  close STDOUT;
  open STDOUT, '>', \$tmp;

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
      'AFR' => '0.234',
      'AMR' => '0.0432',
      'EAS' => '0.0476',
      'EUR' => '0.0467',
      'SAS' => '0.0194',
    },
    'run with haplotype_frequencies'
  );


  # haplotypes no header tests
  $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({%$cfg_hash, output_file => 'STDOUT'});
  $tmp = undef;
  close STDOUT;
  open STDOUT, '>', \$tmp;
  $runner->haplotype_frequencies($test_cfg->{haplo_freqs_nh});
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
      'unknown1' => 'ENST00000307301',
      'freq1' => '0.0905',
      'freq2' => '0.234',
      'freq3' => '0.0432',
      'freq4' => '0.0476',
      'freq5' => '0.0467',
      'freq6' => '0.0194',
    },
    'run with haplotype_frequencies - no header'
  );


  # use input_data
  open IN, $test_cfg->{test_vcf};
  my @lines = <IN>;
  $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({
    %{$test_cfg->base_testing_cfg},
    input_data => join("", @lines),
    output_file => 'STDOUT'
  });
  $tmp = undef;
  close STDOUT;
  open STDOUT, '>', \$tmp;
  $runner->run;

  $tmp =~ s/\t/  /g;
  is(
    $tmp,
  qq{ENST00000352957  ENST00000352957:612A>G    ENSP00000284967:REF      rs1135618  HG00096:1
ENST00000352957  ENST00000352957:91T>C,612A>G    ENSP00000284967:31S>P      rs3989369,rs1135618  HG00096:1
ENST00000307301  ENST00000307301:612A>G    ENSP00000305682:REF      rs1135618  HG00096:1
ENST00000307301  ENST00000307301:91T>C,612A>G    ENSP00000305682:31S>P      rs3989369,rs1135618  HG00096:1
ENST00000419219  ENST00000419219:582A>G    ENSP00000404426:REF      rs1135618  HG00096:1
ENST00000419219  ENST00000419219:91T>C,582A>G    ENSP00000404426:31S>P      rs3989369,rs1135618  HG00096:1
},
    'run - use input_data'
  );



  # restore STDOUT
  open(STDOUT, ">&SAVE") or die "Can't restore STDOUT\n";
}


## test JSON output
SKIP: {

  use_ok('Bio::EnsEMBL::VEP::Haplo::Runner');

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  no warnings 'once';
  skip 'Set::IntervalTree or JSON not installed', 1 unless $Bio::EnsEMBL::VEP::Haplo::Runner::CAN_USE_JSON;

  my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({%$cfg_hash, output_file => 'STDOUT', json => 1, dont_export => 'seq,aligned_sequences'});

  no warnings 'once';
  open(SAVE, ">&STDOUT") or die "Can't save STDOUT\n"; 

  my $tmp;
  close STDOUT;
  open STDOUT, '>', \$tmp;

  $runner->run();

  my $json = JSON->new();

  is_deeply(
    $json->decode((split("\n", $tmp))[0]),
    {
      'total_population_counts' => {
        '_all' => 2
      },
      'protein_haplotypes' => [
        {
          'count' => 1,
          'flags' => [],
          'population_counts' => {
            '_all' => 1
          },
          'name' => 'ENSP00000284967:REF',
          'other_hexes' => {
            'e5745909c573dfc830946d4e6994a47d' => 1
          },
          'frequency' => '0.5',
          'samples' => {
            'HG00096' => 1
          },
          'contributing_variants' => [],
          'type' => 'protein',
          'population_frequencies' => {
            '_all' => '0.5'
          },
          'hex' => '7aa776c0543c7d064e0d146b60541843',
          'diffs' => [],
          'has_indel' => 0
        },
        {
          'count' => 1,
          'flags' => [],
          'population_counts' => {
            '_all' => 1
          },
          'name' => 'ENSP00000284967:31S>P',
          'other_hexes' => {
            '2ff7fea7ee58c8a776fc97f01f9d5296' => 1
          },
          'frequency' => '0.5',
          'samples' => {
            'HG00096' => 1
          },
          'contributing_variants' => [
            'rs3989369'
          ],
          'type' => 'protein',
          'population_frequencies' => {
            '_all' => '0.5'
          },
          'hex' => 'c6c26cfaa79a2046a684220c04b5b2e4',
          'diffs' => [
            {
              'sift_score' => 1,
              'sift_prediction' => 'tolerated',
              'polyphen_prediction' => 'benign',
              'diff' => '31S>P',
              'polyphen_score' => 0
            }
          ],
          'has_indel' => 0
        }
      ],
      'cds_haplotypes' => [
        {
          'count' => 1,
          'population_counts' => {
            '_all' => 1
          },
          'name' => 'ENST00000352957:91T>C,612A>G',
          'other_hexes' => {
            'c6c26cfaa79a2046a684220c04b5b2e4' => 1
          },
          'frequency' => '0.5',
          'samples' => {
            'HG00096' => 1
          },
          'contributing_variants' => [
            'rs3989369',
            'rs1135618'
          ],
          'type' => 'cds',
          'hex' => '2ff7fea7ee58c8a776fc97f01f9d5296',
          'population_frequencies' => {
            '_all' => '0.5'
          },
          'diffs' => [
            {
              'variation_feature' => 'rs3989369',
              'variation_feature_id' => undef,
              'diff' => '91T>C'
            },
            {
              'diff' => '612A>G'
            }
          ],
          'has_indel' => 0
        },
        {
          'count' => 1,
          'population_counts' => {
            '_all' => 1
          },
          'name' => 'ENST00000352957:612A>G',
          'other_hexes' => {
            '7aa776c0543c7d064e0d146b60541843' => 1
          },
          'frequency' => '0.5',
          'samples' => {
            'HG00096' => 1
          },
          'contributing_variants' => [
            'rs1135618'
          ],
          'type' => 'cds',
          'hex' => 'e5745909c573dfc830946d4e6994a47d',
          'population_frequencies' => {
            '_all' => '0.5'
          },
          'diffs' => [
            {
              'diff' => '612A>G'
            }
          ],
          'has_indel' => 0
        }
      ],
      'transcript_id' => 'ENST00000352957',
      'total_haplotype_count' => 2
    },
    'JSON output'
  );

  # Test empty output (when none of the variants overlap a transcript)
  open IN, $test_cfg->{no_trans_vcf} ;
  my @new_lines = <IN>;
  $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({
    %{$test_cfg->base_testing_cfg},
    input_data => join("", @new_lines),
    output_file => 'STDOUT'
  });
  # Catch and evaluate the warning message
  my $warning = warning { $runner->run() };
  like($warning, qr/Haplosaurus can't find transcripts/, 'warning message because of empty output');
}




done_testing();
