# Copyright [2016] EMBL-European Bioinformatics Institute
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
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation');

my $dir = $test_cfg->{sereal_dir};

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
ok($cfg, 'get new config object');

my $c = Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation->new({
  config => $cfg,
});
ok($c, 'new is defined');


## METHODS
##########

ok($c->filter_variation({failed => 0}),  'filter_variation pass');
ok(!$c->filter_variation({failed => 1}), 'filter_variation fail');

$c->{failed} = 1;
ok($c->filter_variation({failed => 1}), 'filter_variation pass with failed on');
$c->{failed} = 0;


## NOVEL TESTS
##############

use_ok('Bio::EnsEMBL::Variation::VariationFeature');

my $input = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
  start  => 10,
  end    => 10,
  strand => 1,
  allele_string => 'A/G',
});

my $existing = {
  start  => 10,
  end    => 10,
  strand => 1,
  allele_string => 'A/G',
};

ok(!$c->is_var_novel($existing, $input), 'is_var_novel exact match');

$existing->{allele_string} = 'A/T';
ok($c->is_var_novel($existing, $input), 'is_var_novel alleles dont match but no check');

$c->{no_check_alleles} = 1;
ok(!$c->is_var_novel($existing, $input), 'is_var_novel alleles dont match with check');

$existing = {
  start  => 10,
  end    => 10,
  strand => -1,
  allele_string => 'T/C',
};
ok(!$c->is_var_novel($existing, $input), 'is_var_novel rev strand exact match');

$c->{no_check_alleles} = 0;
$existing->{allele_string} = 'T/G';
ok($c->is_var_novel($existing, $input), 'is_var_novel rev strand alleles dont match but no check');

$c->{no_check_alleles} = 1;
ok(!$c->is_var_novel($existing, $input), 'is_var_novel rev strand alleles dont match with check');



# done
done_testing();
