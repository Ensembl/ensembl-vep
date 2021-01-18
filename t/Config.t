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
ok(my $tmp = $cfg->read_config_from_file($test_cfg->{test_ini_file}), 'read_config_from_file ok');

is($tmp->{test1}, 'hello', 'read_config_from_file basic');
is($tmp->{test2}, 'foo bar', 'quoted string with spaces');
is($tmp->{test3}, 'foo', 'read_config_from_file flag not allowed multiple');
is($tmp->{individual}, 'dave,barry', 'read_config_from_file flag list preserved comma-separated');
is_deeply($tmp->{plugin}, [qw(foo bar too)], 'read_config_from_file flag allowed multiple');

$cfg->read_config_from_file($test_cfg->{test_ini_file}, $tmp);
is_deeply($tmp->{plugin}, [qw(foo bar too foo bar too)], 'read_config_from_file flag multiples added not overwritten');





## new() METHOD LOGIC
#####################

# defaults applied
$cfg = Bio::EnsEMBL::VEP::Config->new();
is($cfg->param('species'), 'homo_sapiens', 'defaults applied');

# empty arrayref type does not trigger option set
$cfg = Bio::EnsEMBL::VEP::Config->new();
is($cfg->param('regulatory'), undef, 'empty array ref (cell_type) does not trigger option set');

# option sets
$cfg = Bio::EnsEMBL::VEP::Config->new({genomes => 1});
ok(($cfg->param('host') eq 'mysql-eg-publicsql.ebi.ac.uk' and $cfg->param('port') == 4157), 'option sets, one in multiple out');

$cfg = Bio::EnsEMBL::VEP::Config->new({af => 1});
is($cfg->param('check_existing'), 1, 'option sets, multiple in same out 1');

$cfg = Bio::EnsEMBL::VEP::Config->new({check_frequency => 1});
is($cfg->param('check_existing'), 1, 'option sets, multiple in same out 2');

$cfg = Bio::EnsEMBL::VEP::Config->new({gff => 'test'});
is_deeply($cfg->param('custom'), ['test,,gff'], 'option sets, substitution');

$cfg = Bio::EnsEMBL::VEP::Config->new({ucsc_assembly => 'hg38', phyloP => [7, 100], custom => []});
is_deeply(
  $cfg->param('custom'),
  [
    'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP7way/hg38.phyloP7way.bw,phlyoP7way,bigwig,exact',
    'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw,phlyoP100way,bigwig,exact'
  ],
  'option sets, multiple substitution'
);

# give config file
$cfg = Bio::EnsEMBL::VEP::Config->new({config => $test_cfg->{test_ini_file}});
is($cfg->param('test1'), 'hello', 'config file');

# give config file auto-detected as $config->{dir}.'/vep.ini'
$cfg = Bio::EnsEMBL::VEP::Config->new({dir => $test_cfg->{cache_root_dir}.'/../'});
is($cfg->param('test1'), 'hello', 'ini file');

# list conversion
$cfg = Bio::EnsEMBL::VEP::Config->new({individual => 'dave,barry,keith'});
is_deeply($cfg->param('individual'), [qw(dave barry keith)], 'list conversion');

# deprecated
throws_ok { Bio::EnsEMBL::VEP::Config->new({convert => 1}) } qr/deprecated/, 'deprecated no replacement';
#throws_ok { Bio::EnsEMBL::VEP::Config->new({gmaf => 1}) } qr/deprecated.+\-\-af/, 'deprecated with replacement';

# invalid value
throws_ok { Bio::EnsEMBL::VEP::Config->new({format => 'gobbledegook'}) } qr/not a valid value/, 'invalid value for key';

# incompatible flags
throws_ok { Bio::EnsEMBL::VEP::Config->new({database => 1, cache => 1}) } qr/Can\'t use.+together/, 'incompatible params';

# incompatible ok with safe on
ok(Bio::EnsEMBL::VEP::Config->new({safe => 1, most_severe => 1, symbol => 1}), 'incompatible pass with safe');

# required
throws_ok { Bio::EnsEMBL::VEP::Config->new({phyloP => 1}) } qr/You must set --\w+ to use --\w+/, 'required params';

# missing database/cache/offline
throws_ok { Bio::EnsEMBL::VEP::Config->new({database => 0}) } qr/The VEP can read gene data from/, 'no database/cache/offline';


# some specific case tests
$cfg = Bio::EnsEMBL::VEP::Config->new({output_file => 'STDOUT', verbose => 1});
is($cfg->param('quiet'), 1, 'STDOUT output turns on quiet');
is($cfg->param('verbose'), undef, 'STDOUT output turns off verbose');

$cfg = Bio::EnsEMBL::VEP::Config->new({everything => 1, database => 1});
is($cfg->param($_), undef, 'everything with database turns off '.$_) for qw(af_1kg af_esp af_exac af_gnomad pubmed);

#$cfg = Bio::EnsEMBL::VEP::Config->new({species => 'mus_musculus1'});
#is($cfg->param('species'), 'mus_musculus', '1 getting trimmed from species name');



## DONE
#######
done_testing();

