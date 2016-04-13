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
is_deeply($ib->_full_buffer, [], '_full_buffer');

push @{$ib->buffer}, 'hello';
$ib->reset_buffer;
is_deeply($ib->buffer, [], 'reset_buffer');

$vfs = $ib->next();
is(ref($vfs), 'ARRAY', 'next ref');
is(scalar @$vfs, $ib->param('buffer_size'), 'next size');

delete $vfs->[0]->{adaptor};
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

delete $vfs->[0]->{adaptor};
is_deeply($vfs->[0], bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'rs148490508',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => 25592911,
  'start' => '25592911'
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'next again first variant');


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


$vfs = $ib->next();
is(scalar @$vfs, 0, 'next again - finished');

$ib->reset_buffer();
delete($ib->{_config});
is_deeply($ib, bless( {
  'buffer_size' => 5,
  '_full_buffer' => [],
  '_buffer' => []
}, 'Bio::EnsEMBL::VEP::InputBuffer' ), 'finished buffer empty after reset_buffer');



# done
done_testing();
