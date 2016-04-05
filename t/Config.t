# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new();
ok($cfg, 'new is defined');

is(ref($cfg), 'Bio::EnsEMBL::VEP::Config', 'check class');

throws_ok { Bio::EnsEMBL::VEP::Config->new('not a hashref') } qr/First argument/, 'new with invalid arg';



## METHODS
##########

$cfg = Bio::EnsEMBL::VEP::Config->new({testing => 'hello'});
is($cfg->param('testing'), 'hello', 'param get');

is($cfg->param('testing', 'goodbye'), 'goodbye', 'param set');

throws_ok { $cfg->param() } qr/No param/, 'param without key';

# ini reading
$cfg = Bio::EnsEMBL::VEP::Config->new();

throws_ok { $cfg->read_config_from_file() } qr/Could not open config file/, 'read_config_from_file no file given';
throws_ok { $cfg->read_config_from_file('does_not_exist') } qr/Could not open config file/, 'read_config_from_file invalid file given';
ok(my $tmp = $cfg->read_config_from_file($Bin.'/testdata/test_vep.ini'), 'read_config_from_file ok');

is($tmp->{test1}, 'hello', 'read_config_from_file basic');
is($tmp->{test2}, 'foo bar', 'quoted string with spaces');
is($tmp->{test3}, 'foo', 'read_config_from_file flag not allowed multiple');
is($tmp->{individual}, 'dave,barry', 'read_config_from_file flag list preserved comma-separated');
is_deeply($tmp->{plugin}, [qw(foo bar too)], 'read_config_from_file flag allowed multiple');





## new() METHOD LOGIC
#####################

# defaults applied
$cfg = Bio::EnsEMBL::VEP::Config->new();
is($cfg->param('species'), 'homo_sapiens', 'defaults applied');

# option sets
$cfg = Bio::EnsEMBL::VEP::Config->new({genomes => 1});
ok(($cfg->param('host') eq 'mysql-eg-publicsql.ebi.ac.uk' and $cfg->param('port') == 4157), 'option sets, one in multiple out');

$cfg = Bio::EnsEMBL::VEP::Config->new({gmaf => 1});
is($cfg->param('check_existing'), 1, 'option sets, multiple in same out 1');

$cfg = Bio::EnsEMBL::VEP::Config->new({check_alleles => 1});
is($cfg->param('check_existing'), 1, 'option sets, multiple in same out 2');

# list conversion
$cfg = Bio::EnsEMBL::VEP::Config->new({individual => 'dave,barry,keith'});
is_deeply($cfg->param('individual'), [qw(dave barry keith)], 'list conversion');

# invalid value
throws_ok { Bio::EnsEMBL::VEP::Config->new({format => 'gobbledegook'}) } qr/not a valid value/, 'invalid value for key';

# incompatible flags
throws_ok { Bio::EnsEMBL::VEP::Config->new({database => 1, cache => 1}) } qr/Can\'t use.+together/, 'incompatible params';

# some specific case tests
$cfg = Bio::EnsEMBL::VEP::Config->new({output_file => 'STDOUT', verbose => 1});
is($cfg->param('quiet'), 1, 'STDOUT output turns on quiet');
is($cfg->param('verbose'), undef, 'STDOUT output turns off verbose');

$cfg = Bio::EnsEMBL::VEP::Config->new({everything => 1, database => 1});
is($cfg->param($_), undef, 'everything with database turns off '.$_) for qw(maf_1kg maf_esp maf_exac pubmed);

throws_ok { Bio::EnsEMBL::VEP::Config->new({tabix => 1}) } qr/Output must be vcf/, 'output must be vcf to use --tabix';

throws_ok { Bio::EnsEMBL::VEP::Config->new({original => 1}) } qr/provide output filters/, 'must provide --filters with --original';



## DONE
#######
done_testing();

