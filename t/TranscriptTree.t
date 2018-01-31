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
  skip 'Set::IntervalTree not installed', 37 unless $Bio::EnsEMBL::VEP::AnnotationType::Transcript::CAN_USE_INTERVAL_TREE;

  # use test
  use_ok('Bio::EnsEMBL::VEP::TranscriptTree');

  throws_ok {Bio::EnsEMBL::VEP::TranscriptTree->new()} qr/reference.+undef/, 'new - no annotation source';

  throws_ok {
    Bio::EnsEMBL::VEP::TranscriptTree->new({annotation_source => bless({}, 'test')})
  } qr/not an ISA of .+Transcript/, 'new - wrong class';

  use_ok('Bio::EnsEMBL::VEP::Haplo::Runner');
  ok(my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new($test_cfg->base_testing_cfg), 'get runner');

  my $t;
  ok($t = Bio::EnsEMBL::VEP::TranscriptTree->new({annotation_source => $runner->get_all_AnnotationSources->[0], config => $runner->config}), 'new');

  is(ref($t), 'Bio::EnsEMBL::VEP::TranscriptTree', 'check class');


  ## METHOD TESTS
  ###############

  is_deeply($t->valid_chromosomes, [21, 'LRG_485'], 'valid_chromosomes - get');
  is_deeply($t->valid_chromosomes([21]), [21], 'valid_chromosomes - set');

  is(ref($t->get_chr_tree('foo')), 'Set::IntervalTree', 'get_chr_tree - ref');

  $t->insert('foo', 1, 5);
  is_deeply($t->fetch('foo', 2, 3), [[1, 5]], 'fetch 1');
  is_deeply($t->fetch('foo', 0, 1), [[1, 5]], 'fetch 2');
  is_deeply($t->fetch('foo', 0, 6), [[1, 5]], 'fetch 3');
  is_deeply($t->fetch('foo', 4, 6), [[1, 5]], 'fetch 4');
  is_deeply($t->fetch('foo', 6, 7), [], 'fetch 5');


  $t->insert('foo', 4, 8);
  is_deeply($t->fetch('foo', 2, 3), [[1, 8]], 'fetch after update 1');

  $t->insert('foo', 9, 10);
  is_deeply($t->fetch('foo', 2, 9), [[1, 8], [9, 10]], 'fetch after update 2');

  ok($t->chromosome_synonyms($test_cfg->{chr_synonyms}), 'load synonyms');

  is_deeply($t->fetch('NC_000021.9', 25585733, 25585733), [[25585656, 25607517]], 'fetch - use synonyms');

  is_deeply($t->fetch('chr21', 25585733, 25585733), [[25585656, 25607517]], 'fetch - remove chr');

  $t = Bio::EnsEMBL::VEP::TranscriptTree->new({annotation_source => $runner->get_all_AnnotationSources->[0], config => $runner->config});
  $t->valid_chromosomes(['chr21']);
  $t->insert('chr21', 5, 10);
  is_deeply($t->fetch(21, 4, 6), [[5, 10]], 'fetch - add chr');

  # insert object
  $t->insert('chrobj', 5, 10, {foo => 'bar'});
  is_deeply($t->fetch('chrobj', 7, 8), [{foo => 'bar', s => 5, e => 10}], 'insert and fetch object');

  # insert another doesn't merge
  $t->insert('chrobj', 6, 12, {goo => 'car'});
  is_deeply($t->fetch('chrobj', 7, 8), [{foo => 'bar', s => 5, e => 10}, {goo => 'car', s => 6, e => 12}], 'insert another doesn\'t merge');


  ## _get_obj_start_end util method
  is_deeply($t->_get_obj_start_end([1, 5]), [1, 5], '_get_obj_start_end - arrayref');
  is_deeply($t->_get_obj_start_end({s => 1, e => 5}), [1, 5], '_get_obj_start_end - hashref 1');
  is_deeply($t->_get_obj_start_end({start => 1, end => 5}), [1, 5], '_get_obj_start_end - hashref 2');
  is_deeply($t->_get_obj_start_end({s => 1}), [1, 1], '_get_obj_start_end - no end');
  throws_ok {$t->_get_obj_start_end({})} qr/No start field/, '_get_obj_start_end throws no start';

  ## _get_dist util method
  is($t->_get_dist(1, 5, 6, 7), 1, '_get_dist arrayref 1');
  is($t->_get_dist(6, 8, 2, 4), 2, '_get_dist arrayref 2');
  is($t->_get_dist(1, 5, 4, 5), 0, '_get_dist arrayref overlap 0');
  is($t->_get_dist(1, 5, -3, -1), 2, '_get_dist arrayref negative');


  ## nearest
  is_deeply($t->nearest('nearest', 1, 5), [], 'nearest - no result');

  $t->insert('nearest', 1, 5, {foo => 1});
  $t->insert('nearest', 11, 16, {bar => 1});
  is_deeply($t->nearest('nearest', 6, 6), [{foo => 1, s => 1, e => 5}], 'nearest - left');
  is_deeply($t->nearest('nearest', 9, 9), [{bar => 1, s => 11, e => 16}], 'nearest - right');
  is_deeply($t->nearest('nearest', 8, 8), [{foo => 1, s => 1, e => 5}, {bar => 1, s => 11, e => 16}], 'nearest - 2 hits');
  is_deeply($t->nearest('nearest', 7, 9), [{foo => 1, s => 1, e => 5}, {bar => 1, s => 11, e => 16}], 'nearest - 2 hits input not length 1');
}


done_testing();
