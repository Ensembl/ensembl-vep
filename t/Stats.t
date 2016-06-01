# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
    }
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
    }
  },
  'log_VariationFeature 2'
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
