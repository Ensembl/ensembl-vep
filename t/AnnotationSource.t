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
my $test_cfg = VEPTestingConfig->new();

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource');

my $dir = $test_cfg->{cache_dir};

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
ok($cfg, 'get new config object');

my $as = Bio::EnsEMBL::VEP::AnnotationSource->new({
  config => $cfg,
  cache_region_size => 1000000,
  up_down_size => 0,
  info => { foo => 'bar' }
});
ok($as, 'new is defined');

is_deeply($as->info, { foo => 'bar' }, 'info');



## TESTS WITH AN INPUT BUFFER
#############################

use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21, 1]});
ok($p, 'get parser object');

use_ok('Bio::EnsEMBL::VEP::InputBuffer');
my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'check class');

is(ref($ib->next()), 'ARRAY', 'check buffer next');

is_deeply(
  $as->get_all_regions_by_InputBuffer($ib),
  [[21, 25]],
  'get_all_regions_by_InputBuffer'
);

is_deeply(
  $ib->min_max,
  [25585733, 25982445],
  'InputBuffer min_max'
);


my $inputs = [
  [ [qw(1 999999 . A G . . .)],    [[1, 0]],         'before boundary 1'],
  [ [qw(1 1000000 . A G . . .)],   [[1, 0]],         'before boundary 2'],
  [ [qw(1 1000001 . A G . . .)],   [[1, 1]],         'after boundary'],
  [ [qw(1 1000000 . AC GT . . .)], [[1, 0], [1, 1]], 'across boundary'],
  [ [qw(1 1000000 . C CT . . .)],  [[1, 0], [1, 1]], 'across boundary insertion'],

  [
    [qw(1 1 . . <DEL> . . SVTYPE=DEL;END=3000001)],
    [[1, 0], [1, 1], [1, 2], [1, 3]],
    'spans mulitple'
  ],
  [
    [[qw(1 1000000 . A G . . .)], [qw(1 1000001 . A G . . .)]],
    [[1, 0], [1, 1]],
    'one before, one after'
  ],
];

is_deeply(get_regions_from_input($_->[0]), $_->[1], $_->[2] || join(" ", @{$_->[0]})) for @$inputs;

# done
done_testing();



sub get_regions_from_input {
  my $input = shift;

  my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);

  my $as = Bio::EnsEMBL::VEP::AnnotationSource->new({
    config => $cfg,
    cache_region_size => 1000000,
    up_down_size => 0,
  });

  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      file => $test_cfg->create_input_file($input),
      valid_chromosomes => [1]
    })
  });
  $ib->next();

  return $as->get_all_regions_by_InputBuffer($ib);
}
