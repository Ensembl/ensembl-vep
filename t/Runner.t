# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
use B;

use lib $Bin;
use VEPTestingConfig;

my $CAN_USE_CAPTURE_TINY;
BEGIN {
  if(eval qq{require Capture::Tiny; use Capture::Tiny qw(capture); 1}) {
    $CAN_USE_CAPTURE_TINY = 1;
  }
}

my $test_cfg = VEPTestingConfig->new();

my $cfg_hash = $test_cfg->base_testing_cfg;
$cfg_hash->{input_data} = '21 25585733 25585733 C/T + rs142513484';

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Runner');

my $runner = Bio::EnsEMBL::VEP::Runner->new($cfg_hash);
ok($runner, 'new is defined');

is(ref($runner), 'Bio::EnsEMBL::VEP::Runner', 'check class');



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
      'gencode_primary' => undef,
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
        'gnomADe' => '170228',
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
    }, 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript' )
  ],
  'get_all_AnnotationSources'
);

# setup_db_connection should return silently in offline mode
ok(!$runner->setup_db_connection(), 'setup_db_connection');

is_deeply($runner->valid_chromosomes, [21, 22, 'LRG_485'], 'valid_chromosomes');

is_deeply($runner->get_Parser, bless({
  '_config' => $runner->config,
  'file' => *Bio::EnsEMBL::VEP::Runner::IN,
  'file_bak' => *Bio::EnsEMBL::VEP::Runner::IN,
  'line_number' => 0,
  'check_ref' => undef,
  'lookup_ref' => undef,
  'chr' => undef,
  'dont_skip' => undef,
  'valid_chromosomes' => {21 => 1, 22 => 1, LRG_485 => 1 },
  'minimal' => undef,
  'lrg' => undef,
  'delimiter' => ' ',
  'max_sv_size' => 10000000,
}, 'Bio::EnsEMBL::VEP::Parser::VEP_input' ), 'get_Parser');

is_deeply($runner->get_InputBuffer, bless({
  '_config' => $runner->config,
  'parser' => bless({
    '_config' => $runner->config,
    'file' => *Bio::EnsEMBL::VEP::Runner::IN,
    'file_bak' => *Bio::EnsEMBL::VEP::Runner::IN,
    'line_number' => 0,
    'check_ref' => undef,
    'lookup_ref' => undef,
    'chr' => undef,
    'dont_skip' => undef,
    'valid_chromosomes' => {21 => 1, 22 => 1, LRG_485 => 1 },
    'minimal' => undef,
    'lrg' => undef,
    'delimiter' => ' ',
    'max_sv_size' => 10000000,
  }, 'Bio::EnsEMBL::VEP::Parser::VEP_input' ),
  'buffer_size' => $runner->param('buffer_size'),
  'minimal' => undef,
  'max_not_ordered_variants' => 100,
  'max_not_ordered_variants_distance' => 100
}, 'Bio::EnsEMBL::VEP::InputBuffer' ), 'get_InputBuffer');

my $info = $runner->get_output_header_info;
ok($info->{time}, 'get_output_header_info - time');
ok($info->{api_version} =~ /^\d+$/, 'get_output_header_info - api_version');
ok($info->{vep_version} =~ /^\d+(\.\d+)?$/, 'get_output_header_info - vep_version');
is(ref($info->{input_headers}), 'ARRAY', 'get_output_header_info - input_headers');

is_deeply(
  $info->{version_data}, 
  {
    'sift' => 'sift5.2.2',
    'polyphen' => '2.2.2',
    '1000genomes' => 'phase3',
    'COSMIC' => '80',
    'ESP' => 'V2-SSA137',
    'gnomADe' => '170228',
    'gencode' => 'GENCODE 24',
    'genebuild' => '2014-07',
    'HGMD-PUBLIC' => '20164',
    'regbuild' => '13.0',
    'dbSNP' => '149',
    'ClinVar' => '201704',
    'assembly' => 'GRCh38.p5'
  },
  'get_output_header_info - version_data'
);

is_deeply($runner->get_OutputFactory, bless( {
  '_config' => $runner->config,
  'uniprot' => undef,
  'xref_refseq' => undef,
  'numbers' => undef,
  'polyphen_analysis' => 'humvar',
  'pick_allele' => undef,
  'coding_only' => undef,
  'refseq' => undef,
  'output_format' => 'vep',
  'no_escape' => undef,
  'total_length' => undef,
  'per_gene' => undef,
  'sift' => undef,
  'appris' => undef,
  'canonical' => undef,
  'symbol' => undef,
  'pick_order' => $runner->param('pick_order'),
  'pick' => undef,
  'no_intergenic' => undef,
  'gene_phenotype' => undef,
  'flag_pick_allele' => undef,
  'process_ref_homs' => undef,
  'most_severe' => undef,
  'terms' => 'SO',
  'flag_pick_allele_gene' => undef,
  'flag_pick' => undef,
  'summary' => undef,
  'polyphen' => undef,
  'af' => undef,
  'biotype' => undef,
  'protein' => undef,
  'domains' => undef,
  'pick_allele_gene' => undef,
  'variant_class' => undef,
  'ccds' => undef,
  'hgvsc' => undef,
  'hgvsp' => undef,
  'hgvsg' => undef,
  'hgvsg_use_accession' => undef,
  'hgvsp_use_prediction' => undef,
  'remove_hgvsp_version' => undef,
  'spdi' => undef,
  'merged' => undef,
  'af_1kg' => undef,
  'tsl' => undef,
  'pubmed' => undef,
  'mane' => undef,
  'mane_select' => undef,
  'mane_plus_clinical' => undef,
  'gencode_primary' => undef,
  'flag_gencode_primary' => undef,
  'spdi'  => undef,
  'header_info' => $info,
  'plugins' => [],
  'no_stats' => 1,
  'allele_number' => undef,
  'show_ref_allele' => undef,
  'max_af' => undef,
  'use_transcript_ref' => undef,
  'af_gnomad' => undef,
  'af_gnomade' => undef,
  'af_gnomadg' => undef,
  'transcript_version' => undef,
  'gene_version' => undef,
  'cell_type' => undef,
  'mirna' => undef,
  'ambiguity' => undef,
  'shift_genomic' => undef,
  'shift_3prime' => undef,
  'vcf_string' => undef,
  'clin_sig_allele' => 1,
  'var_synonyms' => undef,
  'ga4gh_vrs' => undef,
}, 'Bio::EnsEMBL::VEP::OutputFactory::VEP_output' ), 'get_OutputFactory');


ok($runner->init, 'init');

is_deeply(
  $runner->_buffer_to_output($runner->get_InputBuffer),
  [
    join("\t", qw(
      rs142513484
      21:25585733
      T
      ENSG00000154719
      ENST00000307301
      Transcript
      3_prime_UTR_variant
      1122
      - - - - -
      IMPACT=MODIFIER;STRAND=-1
    )),
    join("\t", qw(
      rs142513484
      21:25585733
      T
      ENSG00000154719
      ENST00000352957
      Transcript
      missense_variant
      1033
      991
      331
      A/T
      Gca/Aca
      -
      IMPACT=MODERATE;STRAND=-1
    )),
    join("\t", qw(
      rs142513484
      21:25585733
      T
      ENSG00000260583
      ENST00000567517
      Transcript
      upstream_gene_variant
      -
      -
      -
      -
      -
      -
      IMPACT=MODIFIER;DISTANCE=2407;STRAND=-1
    )),
  ],
  '_buffer_to_output'
);

$runner = Bio::EnsEMBL::VEP::Runner->new($cfg_hash);

is(
  $runner->next_output_line, 
  join("\t", qw(
    rs142513484
    21:25585733
    T
    ENSG00000154719
    ENST00000307301
    Transcript
    3_prime_UTR_variant
    1122
    - - - - -
    IMPACT=MODIFIER;STRAND=-1
  )),
  'next_output_line'
);

# output file
$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, output_file => $test_cfg->{user_file}.'.out'});
is(ref($runner->get_output_file_handle), 'FileHandle', 'get_output_file_handle - ref');

$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, output_file => $test_cfg->{user_file}.'.out'});
throws_ok {$runner->get_output_file_handle} qr/Output file .+ already exists/, 'get_output_file_handle - fail on existing';

$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, output_file => 'stdout'});
is($runner->get_output_file_handle, '*main::STDOUT', 'get_output_file_handle - stdout');

unlink($test_cfg->{user_file}.'.out');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  output_file => $test_cfg->{user_file},
  compress_output => 'gzip'
});
$runner->param('compress_output', 'foo'.$$);

throws_ok {$runner->get_output_file_handle} qr/not found in path/, 'get_output_file_handle - compressed - missing binary';

use_ok('Bio::EnsEMBL::VEP::Utils');
SKIP: {
  no warnings 'once';
  skip 'gzip not in path', 4 unless $Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP;
  
  my $compressed_file = $test_cfg->{user_file}.'.out.gz';

  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    output_file => $compressed_file,
    compress_output => 'gzip'
  });

  my $fh = $runner->get_output_file_handle;
  is(ref($fh), 'FileHandle', 'get_output_file_handle - compressed - ref');

  print $fh "Hello world\n";
  close $fh;

  ok(-e $compressed_file, 'get_output_file_handle - compressed - file exists');
  ok(-B $compressed_file, 'get_output_file_handle - compressed - file is binary');

  open IN, "gzip -dc $compressed_file |";
  my @content = <IN>;
  close IN;
  is($content[0], "Hello world\n", 'get_output_file_handle - compressed - file content');

  unlink($compressed_file);
}

SKIP: {
  no warnings 'once';

  skip 'gzip not in path or Capture::Tiny not installed', 2 unless $Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP && $CAN_USE_CAPTURE_TINY;


  my ($stdout, $stderr, @result) = capture {

    $runner = Bio::EnsEMBL::VEP::Runner->new({
      %$cfg_hash,
      output_file => 'stdout',
      compress_output => 'gzip'
    });

    my $fh = $runner->get_output_file_handle;
    is(ref($fh), 'FileHandle', 'get_output_file_handle - compressed stdout - ref');

    print $fh "Hi\n";

    close $fh;
  };

  ok($stdout, 'get_output_file_handle - compressed stdout - wrote OK');
}

## stats file
#############

$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, output_file => $test_cfg->{user_file}.'.out'});
is(ref($runner->get_stats_file_handle('txt')), 'FileHandle', 'get_stats_file_handle - txt ref');
ok(-e $test_cfg->{user_file}.'.out_summary.txt', 'get_stats_file_handle - txt file exists');

delete $runner->{stats_file_handle};
throws_ok {$runner->get_stats_file_handle('txt')} qr/Stats file .+ already exists/, 'get_stats_file_handle - fail on existing';
unlink($test_cfg->{user_file}.'.out_summary.txt');

is(ref($runner->get_stats_file_handle('html')), 'FileHandle', 'get_stats_file_handle - html ref');
ok(-e $test_cfg->{user_file}.'.out_summary.html', 'get_stats_file_handle - html file exists');
unlink($test_cfg->{user_file}.'.out_summary.html');

$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, stats_file => $test_cfg->{user_file}.'.stats', force_overwrite => 1});
$runner->get_stats_file_handle('txt');
ok(-e $test_cfg->{user_file}.'.stats.txt', 'get_stats_file_handle - extension handling 1');
unlink($test_cfg->{user_file}.'.stats.txt');

$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, stats_file => $test_cfg->{user_file}.'.txt', force_overwrite => 1});
unlink($test_cfg->{user_file}.'.txt');
$runner->get_stats_file_handle('txt');
ok(-e $test_cfg->{user_file}.'.txt', 'get_stats_file_handle - extension handling 2');

unlink($test_cfg->{user_file}.'.html');
$runner->get_stats_file_handle('html');
ok(-e $test_cfg->{user_file}.'.html', 'get_stats_file_handle - extension handling 3');

$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, stats_file => $test_cfg->{user_file}.'_stats', force_overwrite => 1});
unlink($test_cfg->{user_file}.'_stats.txt');
$runner->get_stats_file_handle('txt');
ok(-e $test_cfg->{user_file}.'_stats.txt', 'get_stats_file_handle - extension handling 4');

# clean up
unlink($test_cfg->{user_file}.'.txt');
unlink($test_cfg->{user_file}.'.html');
unlink($test_cfg->{user_file}.'_stats.txt');


## skipped variants file handle
###############################

my $skipped_input = qq{21\t2642609\t.\tG\tA
21\t2827694\trs2376870\tCGTGGATGCGGGGAC\tC\tSVTYPE=DEL;END=2827708
21\tpos123\t.\tG\tA
21\t230710048\trs699\tA\tG
chr402\t89828156\trs201106962\tA\tC
21\t23873048\trs867018739\tZ\tY
21\t23873048\trs867018739\tC\tT
21\t804947\trs569478349\tG\tC
21\t39735004\t.\tA\tC
21\t39735057\trs752685201\tG\tT
21\t39735057\trs752685201\tX\tT
21\t39737938\trs753930790\tG\tA
21\t45687010\trs1964272\tG\tA};

my $skipped_variants_file = $test_cfg->{user_file} . '_skipped.txt';
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  skipped_variants_file => $skipped_variants_file,
  input_data => $skipped_input
});
is(ref($runner->get_skipped_variants_file_handle()), 'FileHandle', 'get_skipped_variants_file_handle');
ok(-e $skipped_variants_file, 'get_skipped_variants_file_handle - file created');

# throw error if file already exists
throws_ok {$runner->get_skipped_variants_file_handle()} qr/Skipped variants file .+ already exists/,
          'get_skipped_variants_file_handle - fail if files already exists';

# overwrite existing file
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  format => 'vcf',
  skipped_variants_file => $skipped_variants_file,
  force_overwrite => 1,
  input_data => $skipped_input
});
ok(-e $skipped_variants_file, 'get_skipped_variants_file_handle - force overwrite');

# test run
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n";

my $tmp;
close STDERR;
open STDERR, '>', \$tmp;
ok($runner->run, 'skipped_variants_file - run');

like($tmp, qr/WARNING/, 'skipped_variants_file - warnings');
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

open IN, $skipped_variants_file;
my @skipped_variants_lines = <IN>;
close IN;
is(scalar @skipped_variants_lines, 4, 'skipped_variants_file - count lines');

is_deeply(
  [grep {!/^\#/} @skipped_variants_lines],
  [
    "[VEP skipped the following variants from *Bio::EnsEMBL::VEP::Runner::IN]\n",
    "Line 3    	21 pos123 . G A                         	Invalid start 'pos123' or end '0' coordinate\n",
    "Line 5    	chr402 89828156 rs201106962 A C         	Chromosome chr402 not found in annotation sources or synonyms; chromosome chr402 does not overlap any features\n",
    "Line 6    	21 23873048 rs867018739 Z Y             	Invalid allele string Z/Y or possible parsing error\n"
  ],
  'skipped_variants_file - lines content'
);

# clean up
unlink($skipped_variants_file);


## run method
#############

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  output_file => $test_cfg->{user_file}.'.out',
  stats_file => $test_cfg->{user_file}.'.txt',
  stats_html => 1,
  stats_text => 1,
  no_stats => 0,
});

ok($runner->run, 'run - ok');

open IN, $test_cfg->{user_file}.'.out';
my @tmp_lines = <IN>;
close IN;
is(scalar @tmp_lines, 41, 'run - count lines');

is_deeply(
  [grep {!/^\#/} @tmp_lines],
  [
    "rs142513484\t21:25585733\tT\tENSG00000154719\tENST00000307301\tTranscript\t3_prime_UTR_variant\t1122\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1\n",
    "rs142513484\t21:25585733\tT\tENSG00000154719\tENST00000352957\tTranscript\tmissense_variant\t1033\t991\t331\tA/T\tGca/Aca\t-\tIMPACT=MODERATE;STRAND=-1\n",
    "rs142513484\t21:25585733\tT\tENSG00000260583\tENST00000567517\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;DISTANCE=2407;STRAND=-1\n",
  ],
  'run - lines content'
);

ok(-e $test_cfg->{user_file}.'.txt', 'dump_stats - text exists');


# test setting up/down distance
is(scalar (grep {/stream/} @tmp_lines), 1, 'run - count up/downstream (default)');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  output_file => $test_cfg->{user_file}.'.out',
  stats_file => $test_cfg->{user_file}.'.txt',
  stats_html => 1,
  stats_text => 1,
  distance => 10000,
  force_overwrite => 1,
});

$runner->run;

open IN, $test_cfg->{user_file}.'.out';
@tmp_lines = <IN>;
close IN;
is(scalar (grep {/stream/} @tmp_lines), 2, 'run - count up/downstream (10000)');
is(scalar (grep {/DISTANCE\=7079/} @tmp_lines), 1, 'run - check distance');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  output_file => $test_cfg->{user_file}.'.out',
  stats_file => $test_cfg->{user_file}.'.txt',
  stats_html => 1,
  stats_text => 1,
  distance => '10000,20000',
  force_overwrite => 1,
  no_stats => 0,
});

$runner->run;

open IN, $test_cfg->{user_file}.'.out';
@tmp_lines = <IN>;
close IN;
is(scalar (grep {/stream/} @tmp_lines), 4, 'run - count up/downstream (1000,2000)');

ok(-e $test_cfg->{user_file}.'.html', 'dump_stats - html exists');

unlink($test_cfg->{user_file}.'.txt');
unlink($test_cfg->{user_file}.'.html');
unlink($test_cfg->{user_file}.'.out');


## test restoration of Slice->seq after run
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  fasta => $test_cfg->{fasta},
  });
$runner->init();

# use of B gleaned from http://www.perlmonks.org/?node_id=413556
my $b = B::svref_2object(\&Bio::EnsEMBL::Slice::seq);
is($b->GV->STASH->NAME, 'Bio::EnsEMBL::Variation::Utils::FastaSequence', 'Slice::seq after init');

$runner->finish();
$b = B::svref_2object(\&Bio::EnsEMBL::Slice::seq);
is($b->GV->STASH->NAME, 'Bio::EnsEMBL::Slice', 'Slice::seq after finish');


## check the situation where an annotation source filters out all variants in the current input buffer
## but there are still variants remaining to be processed from the parser

# the second line should get filtered out here
my $in = qq{21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0
21\t25587758\trs116645811\tG\tA\t.\t.\t.\tGT\t0|0
21\t25587701\trs187353664\tT\tC\t.\t.\t.\tGT\t0|0};

# since we're using buffer_size of 1 this will leave an empty buffer
# runner should be aware of this and still continue on with next buffer fill
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  no_stats => 1,
  buffer_size => 1,
  filter_common => 1,
  freq_freq => 0.0012,
  input_data => $in,
  output_file => $test_cfg->{user_file}.'.out',
  pick => 1,
});
$runner->run();

open IN, $test_cfg->{user_file}.'.out';
@tmp_lines = grep {!/^\#/} <IN>;
close IN;
unlink($test_cfg->{user_file}.'.out');

is_deeply([map {(split("\t", $_))[0]} @tmp_lines], [qw(rs142513484 rs187353664)], 'run with filter check pass over empty buffer');

# now check the same again but using the equivalent forked path through the code
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  no_stats => 1,
  buffer_size => 1,
  filter_common => 1,
  freq_freq => 0.0012,
  fork => 2,
  input_data => $in,
  output_file => $test_cfg->{user_file}.'.out',
  pick => 1,
});

$runner->run();

open IN, $test_cfg->{user_file}.'.out';
@tmp_lines = grep {!/^\#/} <IN>;
close IN;
unlink($test_cfg->{user_file}.'.out');

is_deeply([map {(split("\t", $_))[0]} @tmp_lines], [qw(rs142513484 rs187353664)], 'run with filter check pass over empty buffer - forked');


# plugins
$runner = Bio::EnsEMBL::VEP::Runner->new($cfg_hash);
$runner->{plugins} = undef;
ok($runner->get_all_Plugins, 'get_all_Plugins - plugins may be undefined');

$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, plugin => ['TestPlugin'], quiet => 1});

is_deeply(
  $runner->get_all_Plugins,
  [
    bless( {
      'params' => [],
      'variant_feature_types' => [
        'VariationFeature'
      ],
      'feature_types' => [
        'Transcript'
      ],
      'version' => '2.3',
      'config' => $runner->config->{_params},
    }, 'TestPlugin' )
  ],
  'get_all_Plugins'
);

# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

# plugin doesn't compile
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNoCompile'],
  quiet => 1,
});
ok($runner->get_all_Plugins, 'get_all_Plugins - failed to compile');
ok($tmp =~ /Failed to compile plugin/, 'get_all_Plugins - failed to compile message');

# $runner = Bio::EnsEMBL::VEP::Runner->new({
#   %$cfg_hash,
#   warning_file => 'STDERR',
#   plugin => ['TestPluginNoCompile'],
#   quiet => 1,
#   safe => 1
# });
# throws_ok {$runner->get_all_Plugins} qr/Failed to compile plugin/, 'get_all_Plugins - failed to compile safe die';


# plugin new method fails
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNewFails'],
  quiet => 1,
});
ok($runner->get_all_Plugins, 'get_all_Plugins - new fails');
ok($tmp =~ /Failed to instantiate plugin/, 'get_all_Plugins - new fails message');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNewFails'],
  quiet => 1,
  safe => 1
});
throws_ok {$runner->get_all_Plugins} qr/Failed to instantiate plugin/, 'get_all_Plugins - new fails safe die';


# plugin missing required methods
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNoMethods'],
  quiet => 1,
});
ok($runner->get_all_Plugins, 'get_all_Plugins - missing methods');
ok($tmp =~ /required method/, 'get_all_Plugins - missing methods message');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNoMethods'],
  quiet => 1,
  safe => 1
});
throws_ok {$runner->get_all_Plugins} qr/required method/, 'get_all_Plugins - missing methods safe die';


$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, plugin => ['TestPluginRegFeat'], quiet => 1});
is($runner->param('regulatory'), undef, 'get_all_Plugins - switch on regulatory 1');
$runner->get_all_Plugins;
is($runner->param('regulatory'), 1, 'get_all_Plugins - switch on regulatory 2');



# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";



## post_setup_checks
####################

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  dir => $test_cfg->{cache_root_dir}.'/sereal',
  hgvs => 1,
  offline => 1,
});
throws_ok {$runner->post_setup_checks} qr/Cannot generate HGVS coordinates/, 'post_setup_checks - hgvs + offline with no fasta';

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  dir => $test_cfg->{cache_root_dir}.'/sereal',
  check_ref => 1,
  offline => 1,
});
throws_ok {$runner->post_setup_checks} qr/Cannot check reference sequences/, 'post_setup_checks - check_ref + offline with no fasta';

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  dir => $test_cfg->{cache_root_dir}.'/sereal',
  use_transcript_ref => 1,
  offline => 1,
});
throws_ok {$runner->post_setup_checks} qr/Cannot use transcript reference sequences/, 'post_setup_checks - use_transcript_ref + offline with no fasta';

## status_msg tests require we mess with STDOUT
###############################################

# status_msg prints to STDOUT
no warnings 'once';
open(SAVE, ">&STDOUT") or die "Can't save STDOUT\n"; 

close STDOUT;
open STDOUT, '>', \$tmp;

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  dir => $test_cfg->{cache_root_dir}.'/sereal',
  everything => 1,
  offline => 1,
});
is($runner->param('hgvs'), 1, 'post_setup_checks - everything + offline with no fasta disables hgvs - before');
$runner->post_setup_checks();
is($runner->param('hgvs'), 0, 'post_setup_checks - everything + offline with no fasta disables hgvs - after');
ok($tmp =~ /Disabling --hgvs/, 'post_setup_checks - status_msg for above');
$tmp = '';

foreach my $flag(qw(lrg check_sv check_ref hgvs)) {
  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir}.'/sereal',
    $flag => 1,
    cache => 1,
    offline => 0,
    database => 0,
  });
  $runner->post_setup_checks();
  ok($tmp =~ /Database will be accessed when using --$flag/, 'post_setup_checks - info - '.$flag);
  $tmp = '';
}

open(STDOUT, ">&SAVE") or die "Can't restore STDOUT\n";


## CHECK BAM-EDITED REF MATCH PASSES THROUGH
############################################

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  gff => $test_cfg->{bam_edit_gff},
  bam => $test_cfg->{bam_edit_bam},
  input_data => "21 10524654 . C A",
  pick => 1,
  fasta => $test_cfg->{fasta},
  synonyms => $test_cfg->{chr_synonyms},
});

ok($runner->next_output_line, 'bam-edited transcript ref match passes through');


## FORKING
##########

my $input = [
  [qw(21 25585733 rs142513484 C T . . . GT 0|0)],
  [qw(21 25587701 rs187353664 T C . . . GT 0|0)],
  [qw(21 25587758 rs116645811 G A . . . GT 0|0)],
  [qw(21 25588859 rs199510789 C T . . . GT 0|0)],
];

my $exp = [
  'rs142513484 T ENST00000307301',
  'rs142513484 T ENST00000352957',
  'rs142513484 T ENST00000567517',
  'rs187353664 C ENST00000307301',
  'rs187353664 C ENST00000352957',
  'rs187353664 C ENST00000567517',
  'rs116645811 A ENST00000307301',
  'rs116645811 A ENST00000352957',
  'rs116645811 A ENST00000567517',
  'rs199510789 T ENST00000307301',
  'rs199510789 T ENST00000352957',
  'rs199510789 T ENST00000419219'
];

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->create_input_file($input),
  input_data => undef,
});

my @lines;
while(my $line = $runner->next_output_line) {
  my @split = split("\t", $line);
  push @lines, join(" ", $split[0], $split[2], $split[4]);
}

is_deeply(
  \@lines,
  $exp,
  'fork - check no fork'
);

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->create_input_file($input),
  input_data => undef,
  fork => 2,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  my @split = split("\t", $line);
  push @lines, join(" ", $split[0], $split[2], $split[4]);
}

is_deeply(
  \@lines,
  $exp,
  'fork - check fork (2)'
);

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->create_input_file($input),
  input_data => undef,
  fork => 4,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  my @split = split("\t", $line);
  push @lines, join(" ", $split[0], $split[2], $split[4]);
}

is_deeply(
  \@lines,
  $exp,
  'fork - check fork (4)'
);


# check whole input file
# add HGVS to check that FASTA sequence is being retrieved OK
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->{test_vcf},
  input_data => undef,
  fasta => $test_cfg->{fasta},
  hgvs => 1
});

$exp = [];
while(my $line = $runner->next_output_line) {
  push @$exp, $line;
}

# lets check stats aswell
my $exp_stats = $runner->stats->{stats}->{counters};

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->{test_vcf},
  input_data => undef,
  fasta => $test_cfg->{fasta},
  hgvs => 1,
  fork => 2,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  push @lines, $line;
}

is_deeply(
  \@lines,
  $exp,
  'fork - full file check (2)'
);

is_deeply(
  $runner->stats->{stats}->{counters},
  $exp_stats,
  'fork - stats (2)'
);

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->{test_vcf},
  input_data => undef,
  fasta => $test_cfg->{fasta},
  hgvs => 1,
  fork => 4,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  push @lines, $line;
}

is_deeply(
  \@lines,
  $exp,
  'fork - full file check (4)'
);

is_deeply(
  $runner->stats->{stats}->{counters},
  $exp_stats,
  'fork - stats (4)'
);



# check plugin cache data
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->{test_vcf},
  input_data => undef,
  fork => 2,
  quiet => 1,
  plugin => ['TestPluginCache'],
});

while(my $line = $runner->next_output_line) {
  1;
}

ok($runner->get_all_Plugins->[0]->{cache}, 'fork - plugin cache exists');
is(scalar keys %{$runner->get_all_Plugins->[0]->{cache}}, 132, 'fork - plugin cache check');




# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  fork => 2,
});

$runner->param('warning_file', 'STDERR');

$runner->{_test_warning} = 1;

$runner->next_output_line();
ok($tmp =~ 'TEST WARNING', 'fork - test warning');


$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  fork => 1,
});
$runner->param('warning_file', 'STDERR');
$runner->{_test_die} = 1;

throws_ok {$runner->next_output_line} qr/TEST DIE/, 'fork - test die';


$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  fork => 1,
});
$runner->param('warning_file', 'STDERR');
$runner->{_kill_child} = 'KILL';
throws_ok {$runner->next_output_line} qr/Forked process.+ died.+read-through/, 'fork - kill child';


$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  fork => 1,
});
$runner->param('warning_file', 'STDERR');
$runner->{_kill_self} = 'KILL';
throws_ok {$runner->next_output_line} qr/Forked process.+ died.+read-through/, 'fork - child kill self';

# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";


# warnings
is_deeply($runner->warnings, [], 'warnings - empty');

my $test_warning = qq{WARNING: Test 1

-------------------- EXCEPTION --------------------
MSG: Msg 1
STACK S1
STACK S2


ERROR : Test 2

-------------------- EXCEPTION --------------------
MSG: Msg 2
STACK S1
STACK S2
STACK S3
};

$runner->{_warning_string} = $test_warning;

is_deeply(
  $runner->warnings,
  [
    {
      'msg' => 'Test 1: Msg 1',
      'stack' => 'STACK S1
STACK S2
',
      'type' => 'WARNING'
    },
    {
      'msg' => 'Test 2: Msg 2',
      'stack' => 'STACK S1
STACK S2
STACK S3
',
      'type' => 'ERROR'
    }
  ],
  'warnings - with multiple'
);


# Test to check that after variant is filtered out because all individuals are homozygous for reference
# # vep doesn't stop annotating because of empty output array in Runner.pm
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  input_data => undef,
  input_file => $test_cfg->{test_vcf5},
  dir => $test_cfg->{cache_root_dir},
  output_file => $test_cfg->{user_file}.'.out',
  individual => 'all',
});
ok($runner->run, 'run - ok');
open IN, $test_cfg->{user_file}.'.out';
@tmp_lines = <IN>;
close IN;

is(scalar (grep {/^var5_chr22/} @tmp_lines), 5, 'run - individual all');

unlink($test_cfg->{user_file}.'.out');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  input_data => undef,
  input_file => $test_cfg->{test_vcf6},
  dir => $test_cfg->{cache_root_dir},
  output_file => $test_cfg->{user_file}.'.out',
});
ok($runner->run, 'run - ok');
open IN, $test_cfg->{user_file}.'.out';
@tmp_lines = <IN>;
close IN;

is(scalar (grep {/^var5_chr22/} @tmp_lines), 5, 'run - ref and alt are the same');

unlink($test_cfg->{user_file}.'.out');

### Test individual_zyg
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  input_data => undef,
  input_file => $test_cfg->{test_vcf5},
  dir => $test_cfg->{cache_root_dir},
  output_file => $test_cfg->{user_file}.'.out',
  individual_zyg => 'all',
});
ok($runner->run, 'run - ok');
open IN, $test_cfg->{user_file}.'.out';
@tmp_lines = <IN>;
close IN;

is(scalar (grep {/^rs142513484/} @tmp_lines), 3, 'run - individual_zyg all (homozygous)');
is(scalar (grep {/^var5_chr22/} @tmp_lines), 5, 'run - individual_zyg all (heterozygous)');

unlink($test_cfg->{user_file}.'.out');

## run_rest
###########

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'JSON module not available', 1 unless $Bio::EnsEMBL::VEP::OutputFactory::CAN_USE_JSON;

  $runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, format => 'vcf', delimiter => ' '});
  is_deeply(
    $runner->run_rest('21 25585733 rs142513484 C T . . .'),
    [
      {
        'input' => '21 25585733 rs142513484 C T . . .',
        'assembly_name' => 'GRCh38',
        'end' => 25585733,
        'seq_region_name' => '21',
        'strand' => 1,
        'transcript_consequences' => [
          {
            'gene_id' => 'ENSG00000154719',
            'variant_allele' => 'T',
            'cdna_end' => 1122,
            'consequence_terms' => [
              '3_prime_UTR_variant'
            ],
            'strand' => -1,
            'transcript_id' => 'ENST00000307301',
            'cdna_start' => 1122,
            'impact' => 'MODIFIER'
          },
          {
            'gene_id' => 'ENSG00000154719',
            'cds_start' => 991,
            'variant_allele' => 'T',
            'cdna_end' => 1033,
            'protein_start' => 331,
            'codons' => 'Gca/Aca',
            'cds_end' => 991,
            'consequence_terms' => [
              'missense_variant'
            ],
            'protein_end' => 331,
            'amino_acids' => 'A/T',
            'strand' => -1,
            'transcript_id' => 'ENST00000352957',
            'cdna_start' => 1033,
            'impact' => 'MODERATE'
          },
          {
            'gene_id' => 'ENSG00000260583',
            'variant_allele' => 'T',
            'distance' => 2407,
            'consequence_terms' => [
              'upstream_gene_variant'
            ],
            'strand' => -1,
            'transcript_id' => 'ENST00000567517',
            'impact' => 'MODIFIER'
          }
        ],
        'id' => 'rs142513484',
        'most_severe_consequence' => 'missense_variant',
        'allele_string' => 'C/T',
        'start' => 25585733
      }
    ],
    'run_rest'
  );
}



## DATABASE TESTS
#################

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 6 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens') if $can_use_db;
  
  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
  });

  ok($runner->setup_db_connection(), 'db - setup_db_connection');

  # check it is switching species using aliases
  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'human_vep',
  });
  $runner->setup_db_connection;

  is($runner->species, 'homo_vepiens', 'db - species alias');

  $info = $runner->get_output_header_info;

  is($info->{db_host}, $db_cfg->{host}, 'get_output_header_info - db_host');
  ok($info->{db_name} =~ /homo_vepiens_core.+/, 'get_output_header_info - db_name');
  ok($info->{db_version}, 'get_output_header_info - db_version');

  # GA4GH VRS Allele object
  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    fasta => $test_cfg->{fasta},
    input_data => '21 10524654 . C A',
    %$db_cfg,
    database => 1,
    ga4gh_vrs => 1,
    json => 1,
    offline => 0,
    species => 'human_vep',
  });

  $runner->init();
  is_deeply(
     $runner->run_rest(),
     [{
         "allele_string" => "C/A",
         "assembly_name" => "GRCh38",
         "end" => 10524654,
         "id" => ".",
         "input" =>  "21 10524654 . C A",
         "intergenic_consequences" => [
             {
                 "consequence_terms" => [
                     "intergenic_variant"
                 ],
                 "ga4gh_vrs" =>  {
                     "location" => {
                         "interval"  => {
                             "end" =>  10524654,
                             "start" => 10524653,
                             "type" => "SimpleInterval"
                         },
                         "sequence_id" => "refseq:NC_000021.9",
                         "type" => "SequenceLocation"
                     },
                     "state" => {
                         "sequence" => "A",
                         "type" => "SequenceState"
                     },
                     "type" =>  "Allele"
                 },
                 "impact" =>  "MODIFIER",
                 "variant_allele" => "A"
             }
         ],
         "most_severe_consequence" => "intergenic_variant",
         "seq_region_name" =>  "21",
         "start" =>  10524654,
         "strand" => 1
     }],
     'GA4GH VRS Allele',
  );
};



done_testing();
