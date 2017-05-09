# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Test::More;
use Test::Exception;
use FindBin qw($Bin);

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

## BASIC TESTS
##############

use_ok('Bio::EnsEMBL::VEP::AnnotationType::Transcript');

SKIP: {

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  no warnings 'once';
  skip 'Set::IntervalTree not installed', 9 unless $Bio::EnsEMBL::VEP::AnnotationType::Transcript::CAN_USE_INTERVAL_TREE;

  # use test
  use_ok('Bio::EnsEMBL::VEP::Haplo::InputBuffer');

  throws_ok {Bio::EnsEMBL::VEP::Haplo::InputBuffer->new()} qr/reference for attribute .+ undef/, 'new without transcript_tree';

  use_ok('Bio::EnsEMBL::VEP::Haplo::Runner');
  my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({%{$test_cfg->base_testing_cfg}, input_file => $test_cfg->{test_vcf}});

  ok(my $ib = Bio::EnsEMBL::VEP::Haplo::InputBuffer->new({transcript_tree => $runner->get_TranscriptTree}), 'new');

  is_deeply($ib->transcript_tree, $runner->get_TranscriptTree, 'transcript_tree');

  is($ib->get_max_from_tree(21, 25585733, 25585733), 25607517, 'get_max_from_tree');
  is($ib->get_max_from_tree(21, 25607519, 25607519), 0, 'get_max_from_tree missed');

  $ib = $runner->get_InputBuffer;

  is_deeply(
    $ib->next,
    [
      {
        'alleles' => 'T,C',
        'ids' => [
          'rs1135618'
        ],
        'chr' => '21',
        'gts' => {
          'HG00096' => 'C|C'
        },
        'end' => 25597391,
        'start' => 25597391
      },
      {
        'alleles' => 'A,G',
        'ids' => [
          'rs3989369'
        ],
        'chr' => '21',
        'gts' => {
          'HG00096' => 'A|G'
        },
        'end' => 25606638,
        'start' => 25606638
      }
    ],
    'next'
  );

  is_deeply(
    $ib->next,
    [],
    'next again empty'
  );
}

done_testing();