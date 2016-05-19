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
use Bio::EnsEMBL::Variation::VariationFeature;
my $test_cfg = VEPTestingConfig->new();
my $cfg_hash = $test_cfg->base_testing_cfg;


## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Parser');

my $p = Bio::EnsEMBL::VEP::Parser->new({file => $test_cfg->{test_vcf}});
ok($p, 'new is defined');

is(ref($p), 'Bio::EnsEMBL::VEP::Parser', 'check class');

# need config object to use new with format
use_ok('Bio::EnsEMBL::VEP::Config');
my $cfg = Bio::EnsEMBL::VEP::Config->new({fasta => $test_cfg->{fasta}, offline => 1});



## METHOD TESTS
###############

is($p->file, $test_cfg->{test_vcf}, 'file');
is($p->file('stdin'), '*main::STDIN', 'file STDIN');

throws_ok {$p->file('foo')} qr/File.+does not exist/, 'file not exists';

is($p->line_number, 0, 'line_number get');
is($p->line_number(1), 1, 'line_number set');

is_deeply($p->headers, [], 'headers');

$p->file($test_cfg->{test_vcf});
is($p->detect_format, 'vcf', 'detect_format - VCF');


## VALIDATION CHECKS
####################

is($p->validate_svf(), 1, 'validate_svf - not implemented yet');

is($p->validate_vf(get_vf({allele_string => 'G/C'})), 1, 'validate_vf - ok');

$p->{chr} = ['1'];
is($p->validate_vf(get_vf({allele_string => 'G/C'})), 1, 'validate_vf - chr list include');
is($p->validate_vf(get_vf({allele_string => 'G/C', chr => 2})), 0, 'validate_vf - chr list exclude');
delete($p->{chr});

my $vf = get_vf({allele_string => 'G/C', chr => 'chr1'});
$p->validate_vf($vf);
is($vf->{chr}, '1', 'validate_vf - strip "chr"');

$vf = get_vf({allele_string => 'G/C', chr => 'chromosome'});
$p->validate_vf($vf);
is($vf->{chr}, 'chromosome', 'validate_vf - dont strip "chr" if chromosome');

$vf = get_vf({allele_string => 'G/C', chr => 'CHR_1'});
$p->validate_vf($vf);
is($vf->{chr}, 'CHR_1', 'validate_vf - dont strip "chr" if CHR_');

$vf = get_vf({allele_string => 'G/C', chr => 'M'});
$p->validate_vf($vf);
is($vf->{chr}, 'MT', 'validate_vf - convert M to MT');

$vf = get_vf({allele_string => 'g/c'});
$p->validate_vf($vf);
is($vf->{allele_string}, 'G/C', 'validate_vf - uppercase allele_string');


# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
my $tmp;
open STDERR, '>', \$tmp;

$cfg->param('warning_file', 'STDERR');
$p = Bio::EnsEMBL::VEP::Parser->new({file => $test_cfg->{test_vcf}, config => $cfg});

is($p->validate_vf(get_vf({allele_string => 'G/C', start => 'foo'})), 0, 'validate_vf - start is not number 1');
ok($tmp =~ /coordinate invalid/, 'validate_vf - start is not number 2');

is($p->validate_vf(get_vf({allele_string => 'G/C', end => 'foo'})), 0, 'validate_vf - end is not number 1');
ok($tmp =~ /coordinate invalid/, 'validate_vf - end is not number 2');

is($p->validate_vf(get_vf({allele_string => 'G/C', start => 3, end => 1})), 0, 'validate_vf - start > end+1 1');
ok($tmp =~ /start \> end\+1/, 'validate_vf - start > end+1 2');

is($p->validate_vf(get_vf({allele_string => 'Q/R'})), 0, 'validate_vf - invalid allele string 1');
ok($tmp =~ /Invalid allele string/, 'validate_vf - invalid allele string 2');

is($p->validate_vf(get_vf({allele_string => '-/C'})), 0, 'validate_vf - alleles look like an insertion 1');
ok($tmp =~ /Alleles look like an insertion/, 'validate_vf - alleles look like an insertion 2');

is($p->validate_vf(get_vf({allele_string => 'AA/C', start => 1, end => 1})), 0, 'validate_vf - ref allele length doesnt match coords 1');
ok($tmp =~ /Length of reference allele/, 'validate_vf - ref allele length doesnt match coords 1');

# check_ref
$p->{check_ref} = 1;

is($p->validate_vf(get_vf({allele_string => '-/T', start => 2, end => 1})), 1, 'validate_vf - check_ref insertion');

$vf = get_vf({allele_string => 'C/T', chr => 21, start => 25585733, end => 25585733});
is($p->validate_vf($vf), 1, 'validate_vf - check_ref ok');

$vf->{allele_string} = 'G/T';
is($p->validate_vf($vf), 0, 'validate_vf - check_ref fail');
ok($tmp =~ /Specified reference allele.+does not match Ensembl reference allele/, 'validate_vf - check_ref fail msg');

$vf = get_vf({allele_string => 'CTT/TCC', chr => 21, start => 25585733, end => 25585735});
is($p->validate_vf($vf), 1, 'validate_vf - check_ref long ok');

$p->{check_ref} = 0;

# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";



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
$tmp = '21 25587759 25587759 C/A + test';
open IN, '<', \$tmp;
ok($p->file(*IN), 'use string fed to FH as input data');
is($p->detect_format, 'ensembl', 'detect format of string');
is(<IN>, '21 25587759 25587759 C/A + test', 'FH position reset');



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


sub get_vf {
  my $hashref = shift;

  $hashref->{$_} ||= 1 for qw(chr start end strand);

  return Bio::EnsEMBL::Variation::VariationFeature->new_fast($hashref);
}
