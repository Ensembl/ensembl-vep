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

my $vf;

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Parser::VCF');

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new();
ok($cfg, 'get new config object');

my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
ok($p, 'new is defined');

is(ref($p), 'Bio::EnsEMBL::VEP::Parser::VCF', 'check class');



## METHOD TESTS
###############

is(ref($p->parser), 'Bio::EnsEMBL::IO::Parser::VCF4', 'parser');

$vf = $p->next;
delete($vf->{adaptor});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'rs142513484',
  'map_weight' => 1,
  'allele_string' => 'C/T',
  'end' => 25585733,
  'start' => '25585733'
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'next');

is(ref($p->next), 'Bio::EnsEMBL::Variation::VariationFeature', 'next again');



## OTHER TESTS
##############

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
 config => $cfg,
 file => $test_cfg->create_input_file([qw(21 25587759 sv_dup . <DUP> . . SVTYPE=DUP;END=25587769)])
})->next();
delete($vf->{adaptor});

is_deeply($vf, bless( {
  'outer_end' => undef,
  'chr' => '21',
  'inner_end' => undef,
  'outer_start' => undef,
  'end' => 25587769,
  'inner_start' => undef,
  'strand' => 1,
  'class_SO_term' => 'duplication',
  'variation_name' => 'sv_dup',
  'start' => 25587759
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' ), 'StructuralVariationFeature');




# done
done_testing();