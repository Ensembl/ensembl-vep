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

# use test
use_ok('Bio::EnsEMBL::VEP::TranscriptTree');

throws_ok {Bio::EnsEMBL::VEP::TranscriptTree->new()} qr/reference.+undef/, 'new - no annotation source';

throws_ok {
  Bio::EnsEMBL::VEP::TranscriptTree->new({annotation_source => bless({}, 'test')})
} qr/not an ISA of .+BaseTranscript/, 'new - wrong class';

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

done_testing();