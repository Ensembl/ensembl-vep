# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

my ($vfs, $tmp, $expected);

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::InputBuffer');

# need to get a config object and parser for further tests
use_ok('Bio::EnsEMBL::VEP::Config');
use_ok('Bio::EnsEMBL::VEP::Parser::VCF');

my $cfg = Bio::EnsEMBL::VEP::Config->new({buffer_size => 10});
ok($cfg, 'get new config object');

my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
ok($p, 'get parser object');

my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'check class');



## METHOD_TESTS
###############

is_deeply($ib->buffer, [], '_buffer');
is_deeply($ib->pre_buffer, [], 'pre_buffer');

push @{$ib->buffer}, 'hello';
$ib->reset_buffer;
is_deeply($ib->buffer, [], 'reset_buffer');

$vfs = $ib->next();
is(ref($vfs), 'ARRAY', 'next ref');
is(scalar @$vfs, $ib->param('buffer_size'), 'next size');

delete $vfs->[0]->{adaptor}; delete $vfs->[0]->{_line};
is_deeply($vfs->[0], bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'rs142513484',
  'map_weight' => 1,
  'allele_string' => 'C/T',
  'end' => 25585733,
  'start' => '25585733'
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'next first variant');


$vfs = $ib->next();
is(scalar @$vfs, $ib->param('buffer_size'), 'next again');
is(scalar @$vfs, $ib->param('buffer_size'), 'next again size');

delete $vfs->[0]->{adaptor}; delete $vfs->[0]->{_line};
is_deeply($vfs->[0], bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'rs148490508',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => 25592911,
  'start' => '25592911'
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'next again first variant');

is_deeply($ib->min_max, [25592911, 25603910], 'min_max');

my %tmp = %{$vfs->[0]};
$tmp{start}++;
my $tmp_ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, variation_features => [\%tmp]});
$tmp_ib->next();
is_deeply($tmp_ib->min_max, [25592911, 25592912], 'min_max with insertion');


SKIP: {

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  no warnings 'once';
  skip 'Set::IntervalTree not installed', 2 unless $Bio::EnsEMBL::VEP::InputBuffer::CAN_USE_INTERVAL_TREE;

  is(ref($ib->interval_tree), 'Set::IntervalTree', 'interval_tree');

  is_deeply(
    $ib->interval_tree->fetch(25592910, 25592911),
    [
      bless( {
        'chr' => '21',
        'strand' => 1,
        'variation_name' => 'rs148490508',
        'map_weight' => 1,
        'allele_string' => 'A/G',
        'end' => 25592911,
        'start' => '25592911'
      }, 'Bio::EnsEMBL::Variation::VariationFeature' )
    ],
    'interval_tree fetch'
  );
}

my $exp = [
  bless( {
    'chr' => '21',
    'strand' => 1,
    'variation_name' => 'rs148490508',
    'map_weight' => 1,
    'allele_string' => 'A/G',
    'end' => 25592911,
    'start' => '25592911'
  }, 'Bio::EnsEMBL::Variation::VariationFeature' )
];

is_deeply(
  $ib->get_overlapping_vfs(25592911, 25592911),
  $exp,
  'get_overlapping_vfs 1'
);

is_deeply(
  $ib->get_overlapping_vfs(25592910, 25592911),
  $exp,
  'get_overlapping_vfs 2'
);

is_deeply(
  $ib->get_overlapping_vfs(25592911, 25592912),
  $exp,
  'get_overlapping_vfs 3'
);

is_deeply(
  $ib->get_overlapping_vfs(25592910, 25592912),
  $exp,
  'get_overlapping_vfs 4'
);

is_deeply(
  $ib->get_overlapping_vfs(25592911, 25592910),
  $exp,
  'get_overlapping_vfs 5'
);

is_deeply(
  $ib->get_overlapping_vfs(25592912, 25592911),
  $exp,
  'get_overlapping_vfs 6'
);

is_deeply(
  $ib->get_overlapping_vfs(25592910, 25592910),
  [],
  'get_overlapping_vfs 7'
);

is_deeply(
  $ib->get_overlapping_vfs(25592912, 25592912),
  [],
  'get_overlapping_vfs 8'
);

SKIP: {

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  no warnings 'once';
  skip 'Set::IntervalTree not installed', 8 unless $Bio::EnsEMBL::VEP::InputBuffer::CAN_USE_INTERVAL_TREE;

  my $orig = $Bio::EnsEMBL::VEP::InputBuffer::CAN_USE_INTERVAL_TREE;
  $Bio::EnsEMBL::VEP::InputBuffer::CAN_USE_INTERVAL_TREE = 0;

  delete $ib->{temp}->{interval_tree};

  is_deeply(
    $ib->get_overlapping_vfs(25592911, 25592911),
    $exp,
    'get_overlapping_vfs no tree 1'
  );

  is_deeply(
    $ib->get_overlapping_vfs(25592910, 25592911),
    $exp,
    'get_overlapping_vfs no tree 2'
  );

  is_deeply(
    $ib->get_overlapping_vfs(25592911, 25592912),
    $exp,
    'get_overlapping_vfs no tree 3'
  );

  is_deeply(
    $ib->get_overlapping_vfs(25592910, 25592912),
    $exp,
    'get_overlapping_vfs no tree 4'
  );

  is_deeply(
    $ib->get_overlapping_vfs(25592911, 25592910),
    $exp,
    'get_overlapping_vfs no tree 5'
  );

  is_deeply(
    $ib->get_overlapping_vfs(25592912, 25592911),
    $exp,
    'get_overlapping_vfs no tree 6'
  );

  is_deeply(
    $ib->get_overlapping_vfs(25592910, 25592910),
    [],
    'get_overlapping_vfs no tree 7'
  );

  is_deeply(
    $ib->get_overlapping_vfs(25592912, 25592912),
    [],
    'get_overlapping_vfs no tree 8'
  );

  $Bio::EnsEMBL::VEP::InputBuffer::CAN_USE_INTERVAL_TREE = $orig;
}

# now use those VFs to create from scratch with VFs instead of a parser
$cfg->param('buffer_size', 5);
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, variation_features => $vfs});

is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'new with vfs - check class');

$vfs = $ib->next();
is(ref($vfs), 'ARRAY', 'new with vfs - next ref');
is(scalar @$vfs, $ib->param('buffer_size'), 'new with vfs - next size');

is_deeply($vfs->[0], bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'rs148490508',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => 25592911,
  'start' => '25592911'
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'new with vfs - next first variant');


$vfs = $ib->next();
is(scalar @$vfs, $ib->param('buffer_size'), 'next again');
is(scalar @$vfs, $ib->param('buffer_size'), 'next again size');

$ib->finish_annotation();
is($vfs->[0]->display_consequence, 'intergenic_variant', 'finish_annotation gives intergenic_variant');


$vfs = $ib->next();
is(scalar @$vfs, 0, 'next again - finished');

$ib->reset_buffer();
delete($ib->{_config});
is_deeply($ib, bless( {
  'buffer_size' => 5,
  'pre_buffer' => [],
  'temp' => {},
}, 'Bio::EnsEMBL::VEP::InputBuffer' ), 'finished buffer empty after reset_buffer');



# done
done_testing();
