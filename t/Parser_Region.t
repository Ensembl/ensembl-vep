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
use_ok('Bio::EnsEMBL::VEP::Parser::Region');

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, fasta => $test_cfg->{fasta}});
ok($cfg, 'get new config object');

my $p = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg, file => $test_cfg->create_input_file(['21:25587759-25587759/A'])
});
ok($p, 'new is defined');

throws_ok {
  Bio::EnsEMBL::VEP::Parser::Region->new({
    config => Bio::EnsEMBL::VEP::Config->new($base_testing_cfg),
    file => $test_cfg->create_input_file(['21:25587759-25587759/A'])
  })
} qr/Cannot use Region format in offline mode/, 'new with --offline no --fasta fails';

throws_ok {
  Bio::EnsEMBL::VEP::Parser::Region->new({
    config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, fasta => $test_cfg->{fasta}, check_ref => 1}),
    file => $test_cfg->create_input_file(['21:25587759-25587759/A'])
  })
} qr/Region format is not compatible with --check_ref/, 'new with --check_ref fails';

is(ref($p), 'Bio::EnsEMBL::VEP::Parser::Region', 'check class');


## FORMAT TESTS
###############

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg, file => $test_cfg->create_input_file(['21:25587759-25587759:1/A']),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{slice});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => '21:25587759-25587759:1/A',
  'map_weight' => 1,
  'allele_string' => 'T/A',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => ['21:25587759-25587759:1/A']
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'basic next test');

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg, file => $test_cfg->create_input_file(['21:25587759-25587759/A']),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{slice});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => '21:25587759-25587759/A',
  'map_weight' => 1,
  'allele_string' => 'T/A',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => ['21:25587759-25587759/A']
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'no strand');

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg, file => $test_cfg->create_input_file(['21:25587759-25587759:-1/G']),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{slice});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '-1',
  'variation_name' => '21:25587759-25587759:-1/G',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => ['21:25587759-25587759:-1/G']
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'negative strand (-1)');

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg, file => $test_cfg->create_input_file(['21:25587759-25587758:1/A']),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{slice});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => '21:25587759-25587758:1/A',
  'map_weight' => 1,
  'allele_string' => '-/A',
  'end' => '25587758',
  'start' => '25587759',
  'seq_region_end' => '25587758',
  'seq_region_start' => '25587759',
  '_line' => ['21:25587759-25587758:1/A']
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'insertion');

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg, file => $test_cfg->create_input_file(['21:25587759-25587759:1/-']),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{slice});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => '21:25587759-25587759:1/-',
  'map_weight' => 1,
  'allele_string' => 'T/-',
  'end' => '25587759',
  'start' => '25587759',
  'seq_region_end' => '25587759',
  'seq_region_start' => '25587759',
  '_line' => ['21:25587759-25587759:1/-']
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'deletion');

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg, file => $test_cfg->create_input_file([qw(21:25587759-25587769/DUP)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{slice});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => '1',
  'variation_name' => '21:25587759-25587769/DUP',
  'class_SO_term' => 'duplication',
  'end' => '25587769',
  'start' => '25587759',
  'seq_region_end' => '25587769',
  'seq_region_start' => '25587759',
  '_line' => [qw(21:25587759-25587769/DUP)]
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' ), 'SV dup');


# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$cfg->param('warning_file', 'STDERR');

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg,
  file => $test_cfg->create_input_file([
    [qw(21:25587769-25587767:1/A)],
    [qw(21:25587759-25587759:1/A)],
  ]),
  valid_chromosomes => [21]
})->next();

is($vf->{start}, 25587759, 'skip VF that fails validation');

$cfg->param('dont_skip', 1);

$vf = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $cfg,
  file => $test_cfg->create_input_file([
    [qw(21:25587769-25587767:1/A)],
    [qw(21:25587759-25587759:1/A)],
  ]),
  valid_chromosomes => [21]
})->next();

is($vf->{start}, 25587769, 'dont skip VF that fails validation with dont_skip');

# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

done_testing();
