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
use FindBin qw($Bin);

use lib $Bin;
use VEPTestingConfig;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::VEP::Utils qw(get_time);
my $test_cfg = VEPTestingConfig->new();

my $cfg_hash = $test_cfg->base_testing_cfg;
$cfg_hash->{no_stats} = 0;
$cfg_hash->{output_file} = $test_cfg->{user_file}.'.out';
$cfg_hash->{force_overwrite} = 1;

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Stats');
use_ok('Bio::EnsEMBL::VEP::Runner');

my $s = Bio::EnsEMBL::VEP::Stats->new;
ok($s, 'new is defined');
is_deeply(
  $s,
  bless({
    stats => {
      counters => {}
    },
  }, 'Bio::EnsEMBL::VEP::Stats'),
  'check object'
);


## METHOD TESTS
###############

$s->{info} = {foo => 'bar'};
is_deeply($s->info, {foo => 'bar'}, 'info');

is($s->start_time, get_time(), 'start_time');
is($s->{stats}->{run_time_start}, time(), 'run_time_start');

sleep(1);
ok($s->run_time =~ /^\d+$/, 'run_time');
is($s->end_time, get_time(), 'end_time');

$s->log_lines_read(10);
is($s->{stats}->{lines_read}, 10, 'log_lines_read');

$s->increment_filtered_variants(5);
is_deeply(
  $s->{stats}->{counters},
  {
    filtered_variants => 5
  },
  'increment_filtered_variants'
);



# we need full setup now
my $runner = get_annotated_buffer_runner({input_file => $test_cfg->{test_vcf}});
my $vf = $runner->get_InputBuffer->buffer->[0];

$s = Bio::EnsEMBL::VEP::Stats->new;
$s->log_VariationFeature($vf);
is_deeply(
  $s->{stats}->{counters},
  {
    'classes' => {
      'SNV' => 1
    },
    'var_cons' => {
      'missense_variant' => 1
    },
    'chr' => {
      '21' => {
        '25000000' => 1
      }
    },
    'allele_changes' => {
      'C/T' => 1
    },
    'var_count' => 1,
  },
  'log_VariationFeature 1'
);

$s->log_VariationFeature($runner->get_InputBuffer->buffer->[1]);
is_deeply(
  $s->{stats}->{counters},
  {
    'classes' => {
      'SNV' => 2
    },
    'var_cons' => {
      'missense_variant' => 2
    },
    'chr' => {
      '21' => {
        '25000000' => 2
      }
    },
    'allele_changes' => {
      'T/C' => 1,
      'C/T' => 1
    },
    'var_count' => 2,
  },
  'log_VariationFeature 2'
);

$runner->get_InputBuffer->buffer->[1]->{allele_string} .= '/G';
$s->log_VariationFeature($runner->get_InputBuffer->buffer->[1]);
is_deeply(
  $s->{stats}->{counters},
  {
    'classes' => {
      'SNV' => 3
    },
    'var_cons' => {
      'missense_variant' => 3
    },
    'chr' => {
      '21' => {
        '25000000' => 3
      }
    },
    'allele_changes' => {
      'T/C' => 2,
      'T/G' => 1,
      'C/T' => 1
    },
    'var_count' => 3,
  },
  'log_VariationFeature 3'
);

$s = Bio::EnsEMBL::VEP::Stats->new;
my $vfoa = $runner->get_OutputFactory->get_all_VariationFeatureOverlapAlleles($vf)->[1];
my $hash = $runner->get_OutputFactory->get_all_output_hashes_by_VariationFeature($vf)->[1];
$s->log_TranscriptVariationAllele($vfoa, $hash);

is_deeply(
  $s->{stats}->{counters},
  {
    'protein_pos' => {
      '9' => 1
    },
    'gene' => {
      'ENSG00000154719' => 1
    }
  },
  'log_TranscriptVariationAllele'
);

$s = Bio::EnsEMBL::VEP::Stats->new;
$s->log_VariationFeatureOverlapAllele($vfoa, $hash);

is_deeply(
  $s->{stats}->{counters},
  {
    'protein_pos' => {
      '9' => 1
    },
    'gene' => {
      'ENSG00000154719' => 1
    },
    'transcript' => {
      'ENST00000352957' => 1
    },
    'consequences' => {
      'missense_variant' => 1
    }
  },
  'log_VariationFeatureOverlapAllele'
);

$runner = get_annotated_buffer_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 sv_dup T . . . SVTYPE=DUP;END=25585735)])
});
$vf = $runner->get_InputBuffer->buffer->[0];
$vfoa = $runner->get_OutputFactory->get_all_StructuralVariationOverlapAlleles($vf)->[1];
$hash = $runner->get_OutputFactory->get_all_output_hashes_by_VariationFeature($vf)->[1];

$s = Bio::EnsEMBL::VEP::Stats->new;
$s->log_TranscriptStructuralVariationAllele($vfoa, $hash);

is_deeply(
  $s->{stats}->{counters},
  {
    'protein_pos' => {
      '9' => 1
    },
    'gene' => {
      'ENSG00000154719' => 1
    }
  },
  'log_TranscriptStructuralVariationAllele'
);

$s = Bio::EnsEMBL::VEP::Stats->new;
$s->log_VariationFeatureOverlapAllele($vfoa, $hash);

is_deeply(
  $s->{stats}->{counters},
  {
    'protein_pos' => {
      '9' => 1
    },
    'gene' => {
      'ENSG00000154719' => 1
    },
    'transcript' => {
      'ENST00000352957' => 1
    },
    'consequences' => {
      'feature_elongation' => 1,
      'coding_sequence_variant' => 1
    }
  },
  'log_VariationFeatureOverlapAllele - SV'
);

# test another type that we don't have a specific method for (intergenic)
$runner = get_annotated_buffer_runner({
  input_file => $test_cfg->create_input_file([qw(21 25832817 . C A . . .)])
});
$vf = $runner->get_InputBuffer->buffer->[0];
$vfoa = $runner->get_OutputFactory->get_all_VariationFeatureOverlapAlleles($vf)->[0];
$hash = $runner->get_OutputFactory->get_all_output_hashes_by_VariationFeature($vf)->[0];

$s = Bio::EnsEMBL::VEP::Stats->new;
$s->log_VariationFeatureOverlapAllele($vfoa, $hash);

is_deeply(
  $s->{stats}->{counters},
  {
    'consequences' => {
      'intergenic_variant' => 1
    }
  },
  'log_VariationFeatureOverlapAllele - intergenic'
);

$s = Bio::EnsEMBL::VEP::Stats->new;
$s->log_sift_polyphen('SIFT', 'deleterious');

is_deeply(
  $s->{stats}->{counters},
  {
    'SIFT' => {
      'deleterious' => 1
    }
  },
  'log_sift_polyphen'
);


$runner = get_annotated_buffer_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 . C T . . .)]),
});
$runner->get_OutputFactory->get_all_lines_by_InputBuffer($runner->get_InputBuffer);
$s = $runner->stats;

my $finished = $s->finished_stats;

is_deeply(
  $finished->{general_stats},
  [
    [
      'Lines of input read',
      1
    ],
    [
      'Variants processed',
      1
    ],
    [
      'Variants filtered out',
      0,
    ],
    [
      'Novel / existing variants',
      '-'
    ],
    [
      'Overlapped genes',
      2
    ],
    [
      'Overlapped transcripts',
      3
    ],
    [
      'Overlapped regulatory features',
      '-'
    ]
  ],
  'finished_stats - general_stats'
);

is_deeply(
  [map {$_->{data}} @{$finished->{charts}}],
  [
    {
      'SNV' => 1
    },
    {
      'missense_variant' => 1
    },
    {
      'missense_variant' => 1,
      'upstream_gene_variant' => 1,
      '3_prime_UTR_variant' => 1
    },
    {
      'missense_variant' => 1
    },
    {
      '21' => 1
    },
    {
      '33' => 0,
      '32' => 0,
      '21' => 0,
      '7' => 0,
      '26' => 0,
      '17' => 0,
      '2' => 0,
      '1' => 0,
      '18' => 0,
      '30' => 0,
      '16' => 0,
      '44' => 0,
      '27' => 0,
      '25' => 1,
      '28' => 0,
      '40' => 0,
      '20' => 0,
      '14' => 0,
      '24' => 0,
      '10' => 0,
      '31' => 0,
      '35' => 0,
      '11' => 0,
      '42' => 0,
      '22' => 0,
      '46' => 0,
      '0' => 0,
      '13' => 0,
      '23' => 0,
      '29' => 0,
      '6' => 0,
      '39' => 0,
      '36' => 0,
      '3' => 0,
      '9' => 0,
      '41' => 0,
      '12' => 0,
      '15' => 0,
      '38' => 0,
      '8' => 0,
      '4' => 0,
      '34' => 0,
      '45' => 0,
      '37' => 0,
      '43' => 0,
      '19' => 0,
      '5' => 0
    },
    {
      '40-50%' => 0,
      '70-80%' => 0,
      '50-60%' => 0,
      '90-100%' => 1,
      '10-20%' => 0,
      '20-30%' => 0,
      '30-40%' => 0,
      '00-10%' => 0,
      '60-70%' => 0,
      '80-90%' => 0
    }
  ],
  'finished_stats - charts data'
);


ok($finished->{run_stats}->[0]->[1] =~ /\d+ \(\d+\)/, 'finished_stats - run_stats - version');
is_deeply(
  [map {$_->[0]} @{$finished->{run_stats}}],
  [
    'VEP version (API)',
    'Annotation sources',
    'Species',
    'Command line options',
    'Start time',
    'End time',
    'Run time',
    'Input file',
    'Output file'
  ],
  'finished_stats - run_stats - remaining headers'
);


my $tmp;
open STATS, '>', \$tmp;
ok($s->dump_text(*STATS), 'dump_text');
ok($tmp =~ /^\[VEP run statistics\]/, 'dump_text - content');
close STATS;


open STATS, '>', \$tmp;
ok($s->dump_html(*STATS), 'dump_html');
ok($tmp =~ /^\<html\>/, 'dump_html - content');
close STATS;

SKIP: {
  my $can_use_lint = eval { require HTML::Lint; 1 };
  
  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'HTML::Lint not installed', 1 unless $can_use_lint;

  my $lint = HTML::Lint->new();

  open STATS, '>', \$tmp;
  $s->dump_html(*STATS);

  $lint->parse($tmp);
  $lint->eof();
  ok(!$lint->errors, 'dump_html - HTML::Lint check no errors') or diag(map {$_->as_string."\n"} $lint->errors);

  close STATS;
}

unlink($test_cfg->{user_file}.'.out_summary.html');
unlink($test_cfg->{user_file}.'.out_summary.txt');


sub get_annotated_buffer_runner {
  my $tmp_cfg = shift;

  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir},
    %$tmp_cfg,
  });

  $runner->init;

  my $ib = $runner->get_InputBuffer;
  $ib->next();
  $_->annotate_InputBuffer($ib) for @{$runner->get_all_AnnotationSources};
  $ib->finish_annotation();

  return $runner;
}

done_testing();
