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
my $base_testing_cfg = $test_cfg->base_testing_cfg;

my ($vf, $tmp, $expected);

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Parser::VEP_input');

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($base_testing_cfg);
ok($cfg, 'get new config object');

my $p = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg, file => $test_cfg->create_input_file([qw(21 25587759 25587759 C/A + test)])
});
ok($p, 'new is defined');

is(ref($p), 'Bio::EnsEMBL::VEP::Parser::VEP_input', 'check class');



## FORMAT TESTS
###############

$vf = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg, file => $test_cfg->create_input_file([qw(21 25587759 25587759 C/A + test)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/A',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => [qw(21 25587759 25587759 C/A + test)]
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'basic next test');

$vf = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg, file => $test_cfg->create_input_file([qw(21 25587759 25587759 C/A - test)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '-1',
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/A',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => [qw(21 25587759 25587759 C/A - test)]
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'negative strand (-)');

$vf = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg, file => $test_cfg->create_input_file([qw(21 25587759 25587759 C/A -1 test)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '-1',
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/A',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => [qw(21 25587759 25587759 C/A -1 test)]
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'negative strand (-1)');

$vf = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg, file => $test_cfg->create_input_file([qw(21 25587759 25587759 C/A)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => '21_25587759_C/A',
  'map_weight' => 1,
  'allele_string' => 'C/A',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => [qw(21 25587759 25587759 C/A)]
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'stubby');

$vf = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg, file => $test_cfg->create_input_file([qw(21 25587759 25587769 DUP + test)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => 'test',
  'class_SO_term' => 'duplication',
  'end' => '25587769',
  'start' => '25587759',
  'seq_region_end' => '25587769',
  'seq_region_start' => '25587759',
  '_line' => [qw(21 25587759 25587769 DUP + test)]
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' ), 'SV dup');


# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$cfg->param('warning_file', 'STDERR');

$vf = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg,
  file => $test_cfg->create_input_file([
    [qw(21 foo bar C/A 1)],
    [qw(21 25587759 25587759 C/A 1)],
  ]),
  valid_chromosomes => [21]
})->next();

is($vf->{start}, 25587759, 'skip VF that fails validation');

$cfg->param('dont_skip', 1);

$vf = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $cfg,
  file => $test_cfg->create_input_file([
    [qw(21 foo bar C/A 1)],
    [qw(21 25587759 25587759 C/A 1)],
  ]),
  valid_chromosomes => [21]
})->next();

is($vf->{start}, 'foo', 'dont skip VF that fails validation with dont_skip');

# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

done_testing();
