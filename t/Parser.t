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

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();
my $cfg_hash = $test_cfg->base_testing_cfg;


## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Parser');

my $p = Bio::EnsEMBL::VEP::Parser->new({file => $test_cfg->{test_vcf}});
ok($p, 'new is defined');

is(ref($p), 'Bio::EnsEMBL::VEP::Parser', 'check class');



## METHOD TESTS
###############

is($p->file, $test_cfg->{test_vcf}, 'file');
is($p->file('stdin'), '*main::STDIN', 'file STDIN');

throws_ok {$p->file('foo')} qr/File.+does not exist/, 'file not exists';

is($p->line_number, 0, 'line_number get');
is($p->line_number(1), 1, 'line_number set');

$p->file($test_cfg->{test_vcf});
is($p->detect_format, 'vcf', 'detect_format - VCF');



## FORMAT TESTS
###############

no warnings "qw";
$p->file($test_cfg->create_input_file([qw(21 25587759 test G C,<DEL:*> . . .')]));
is($p->detect_format, 'vcf', 'detect_format - VCF with *');

$p->file($test_cfg->create_input_file([qw(21 25587758 sv_dup T . . . SVTYPE=DUP;END=25587769')]));
is($p->detect_format, 'vcf', 'detect_format - VCF SV 1');

$p->file($test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . .')]));
is($p->detect_format, 'vcf', 'detect_format - VCF SV 2');

$p->file($test_cfg->create_input_file('rs699'));
is($p->detect_format, 'id', 'detect_format - ID');

$p->file($test_cfg->create_input_file('21:g.25585733C>T'));
is($p->detect_format, 'hgvs', 'detect_format - HGVSg');

$p->file($test_cfg->create_input_file('ENST00000352957.8:c.991G>A'));
is($p->detect_format, 'hgvs', 'detect_format - HGVSc');

$p->file($test_cfg->create_input_file('ENSP00000284967.6:p.Ala331Thr'));
is($p->detect_format, 'hgvs', 'detect_format - HGVSp');

$p->file($test_cfg->create_input_file('21 25587759 25587759 C/A + test'));
is($p->detect_format, 'ensembl', 'detect_format - VEP_input');

$p->file($test_cfg->create_input_file('21 25587759 25587759 C/A/G + test'));
is($p->detect_format, 'ensembl', 'detect_format - VEP_input multiple');

$p->file($test_cfg->create_input_file('21 25587759 25587769 DUP + test'));
is($p->detect_format, 'ensembl', 'detect_format - VEP_input SV');

$p->file($test_cfg->create_input_file([qw(chr1 60 T A)]));
is($p->detect_format, 'pileup', 'detect_format - pileup');

$p->file($test_cfg->create_input_file('21 25587759 25587769'));
is($p->detect_format, undef, 'detect_format - incomplete');



## OTHER TESTS
##############

# new using string as input data
my $tmp = '21 25587759 25587759 C/A + test';
open IN, '<', \$tmp;
ok($p->file(*IN), 'use string fed to FH as input data');
is($p->detect_format, 'ensembl', 'detect format of string');
is(<IN>, '21 25587759 25587759 C/A + test', 'FH position reset');


# need config object to use new with format
use_ok('Bio::EnsEMBL::VEP::Config');
my $cfg = Bio::EnsEMBL::VEP::Config->new();

# new with format configured
$p = Bio::EnsEMBL::VEP::Parser->new({config => $cfg, file => $test_cfg->{test_vcf}, format => 'vcf'});
is(ref($p), 'Bio::EnsEMBL::VEP::Parser::VCF', 'new with explicit format');

# new with invalid format
throws_ok {
  Bio::EnsEMBL::VEP::Parser->new({config => $cfg, file => $test_cfg->{test_vcf}, format => 'foo'})
} qr/Unknown or unsupported format/, 'new with unknown format';

# new with format detection
$p = Bio::EnsEMBL::VEP::Parser->new({config => $cfg, file => $test_cfg->{test_vcf}, format => 'guess'});
is(ref($p), 'Bio::EnsEMBL::VEP::Parser::VCF', 'new with format detection');

# new with unsupported format
throws_ok {
  Bio::EnsEMBL::VEP::Parser->new({
    config => $cfg,
    file => $test_cfg->create_input_file([qw(chr1 60 T A)]),
    format => 'guess'
  })
} qr/Unknown or unsupported format/, 'new with unsupported format';


done_testing();