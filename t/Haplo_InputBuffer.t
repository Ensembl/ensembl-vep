# Copyright [2016-2021] EMBL-European Bioinformatics Institute
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
  skip 'Set::IntervalTree not installed', 13 unless $Bio::EnsEMBL::VEP::AnnotationType::Transcript::CAN_USE_INTERVAL_TREE;

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

  my $vfs = $ib->next;

  is(scalar @$vfs, 34, 'next - count');

  is_deeply(
    $vfs->[0],
    {
      'alleles' => 'C,T',
      'record' => [
        '21', '25585733', 'rs142513484', 'C', 'T', '.', '.', '.', 'GT', '0|0'
      ],
      'ids' => [
        'rs142513484'
      ],
      'chr' => '21',
      'end' => 25585733,
      'start' => '25585733'
    },
    'next - first element'
  );

  is_deeply(
    $vfs->[-1],
    {
      'alleles' => 'A,G',
      'record' => [
        '21', '25607474', 'rs141666402', 'A', 'G', '.', '.', '.', 'GT', '0|0'
      ],
      'ids' => [
        'rs141666402'
      ],
      'chr' => '21',
      'end' => 25607474,
      'start' => 25607474
    },
    'next - last element'
  );

  is_deeply(
    $ib->next->[0],
    {
      'alleles' => 'C,A',
      'record' => [
        '21', '25639842', 'rs8132639', 'C', 'A', '.', '.', '.', 'GT', '0|0'
      ],
      'ids' => [
        'rs8132639'
      ],
      'chr' => '21',
      'end' => 25639842,
      'start' => 25639842
    },
    'next again - first element'
  );

  ok(@{$ib->next} > 0, 'next again 2');
  

  is_deeply(
    $ib->next,
    [],
    'next again 3 empty'
  );


  # Test that the input buffer checks the vf chromosome
  my $runner2 = Bio::EnsEMBL::VEP::Haplo::Runner->new({%{$test_cfg->base_testing_cfg}, input_file => $test_cfg->{test_vcf3}});

  $ib = $runner2->get_InputBuffer;

  # Chr 21
  $vfs = $ib->next;
  ok(scalar(@$vfs) == 3 && $vfs->[0]->{chr} == 21, 'The 3 entries from chromosome 21 are in the same buffer');

  # Chr 22
  $vfs = $ib->next;
  ok(scalar(@$vfs) == 1 && $vfs->[0]->{chr} == 22, 'The entry from chromosome 22 is in a different buffer');

  # Test non numeric chromosomes 
  my $runner3 = Bio::EnsEMBL::VEP::Haplo::Runner->new({%{$test_cfg->base_testing_cfg}, input_file => $test_cfg->{test_vcf4}});

  $ib = $runner3->get_InputBuffer;

  # chr21 
  $vfs = $ib->next;
  ok(scalar(@$vfs) == 3 && $vfs->[0]->{chr} eq 'chr21', 'The 3 entries from chr21 are in the same buffer');
  
  # chr22
  $vfs = $ib->next;
  ok(scalar(@$vfs) == 1 && $vfs->[0]->{chr} eq 'chr22', 'The entry from chr22 is in a different buffer');
  
}

done_testing();
