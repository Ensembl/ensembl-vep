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

$p = Bio::EnsEMBL::VEP::Parser->new({file => $test_cfg->{test_vcf}, config => $cfg});
is($p->detect_format, 'vcf', 'detect_format - VCF');

ok($p->file($test_cfg->{test_gzvcf}) =~ /GLOB/, 'gzipped VCF');
is($p->detect_format, 'vcf', 'detect_format - gzipped VCF');

$p->file($test_cfg->{test_vcf});

is($p->delimiter("\t"), "\t", 'delimiter - set');
is($p->delimiter(), "\t", 'delimiter - get');

is_deeply($p->valid_chromosomes, {}, 'valid_chromosomes empty');

$p = Bio::EnsEMBL::VEP::Parser->new({file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
is_deeply($p->valid_chromosomes, {21 => 1}, 'valid_chromosomes');

# allele minimising
is_deeply(
  $p->minimise_alleles([get_vf({allele_string => 'GA/GT'})])->[0],
  bless( {
    'chr' => 1,
    'minimised' => 1,
    'original_allele_string' => 'GA/GT',
    'original_end' => 1,
    'end' => 1,
    'original_start' => 1,
    'strand' => 1,
    'allele_string' => 'A/T',
    'start' => 2,
    'seq_region_start' => 2,
    'seq_region_end' => 1,
  }, 'Bio::EnsEMBL::Variation::VariationFeature'),
  "minimal"
);

is(
  $p->minimise_alleles([get_vf({allele_string => 'GAC/GTC'})])->[0]->{allele_string},
  'A/T',
  "minimal - trim both ends"
);

is_deeply(
  $p->minimise_alleles([get_vf({allele_string => 'TTCCTTCCGACGGTACACACACACA/TTCCTTCCGTCGGTACACACACACA'})])->[0]->{allele_string},
  'A/T',
  "minimal - long trim"
);


## VALIDATION CHECKS
####################

$p = Bio::EnsEMBL::VEP::Parser->new({file => $test_cfg->{test_vcf}, config => $cfg});
$p->{valid_chromosomes} = {1 => 1, 21 => 1, CHR_1 => 1, MT => 1, chromosome => 1, chr12 => 1};

is($p->validate_svf(), 1, 'validate_svf - not implemented yet');

is($p->validate_vf(get_vf({allele_string => 'G/C'})), 1, 'validate_vf - ok');

my $vf = get_vf({allele_string => 'G/C'});
$p->validate_vf($vf);
is($vf->{variation_name}, '1_1_G/C', 'validate_vf - create variation_name');

$p->{chr} = ['1'];
is($p->validate_vf(get_vf({allele_string => 'G/C'})), 1, 'validate_vf - chr list include');
is($p->validate_vf(get_vf({allele_string => 'G/C', chr => 2})), 0, 'validate_vf - chr list exclude');
delete($p->{chr});

$vf = get_vf({allele_string => 'g/c'});
$p->validate_vf($vf);
is($vf->{allele_string}, 'G/C', 'validate_vf - uppercase allele_string');

## have_chr checks
is($p->_have_chr(get_vf({chr => 1})), 1, 'have_chr - pass');
is($p->_have_chr(get_vf({chr => 'chr1'})), 1, 'have_chr - remove chr');
is($p->_have_chr(get_vf({chr => 12})), 1, 'have_chr - add chr');
is($p->_have_chr(get_vf({chr => 2})), 0, 'have_chr - fail');

$vf = get_vf({allele_string => 'G/C', chr => 'M'});
is($p->_have_chr($vf), 1, 'have_chr - M');
is($vf->{chr}, 'MT', 'have_chr - convert M to MT');

ok($p->chromosome_synonyms($test_cfg->{chr_synonyms}), "load chr_synonyms");
$vf = get_vf({chr => 'NC_000021.9'});
is($p->_have_chr($vf), 1, 'have_chr - synonym');
is($vf->{chr}, 'NC_000021.9', 'chr unchanged after synonym check');



# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
my $tmp;
open STDERR, '>', \$tmp;

$cfg->param('warning_file', 'STDERR');
$p = Bio::EnsEMBL::VEP::Parser->new({file => $test_cfg->{test_vcf}, config => $cfg});
$p->{valid_chromosomes} = {1 => 1, 21 => 1};

is($p->validate_vf(get_vf({allele_string => 'A/C', chr => 2})), 0, 'validate_vf - chromosome not in valid list 1');
ok($tmp =~ /Chromosome .* not found/, 'validate_vf - chromosome not in valid list 2');

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

# lookup_ref
$p->{lookup_ref} = 1;
$p->param('fasta', $test_cfg->{fasta});
$p->fasta_db();

$vf = get_vf({chr => 21, start => 25585733, end => 25585733, allele_string => 'N/T'});
is($p->validate_vf($vf), 1, 'validate_vf - lookup_ref pass');
is($vf->allele_string, 'C/T', 'validate_vf - lookup_ref seq OK');

$vf = get_vf({chr => 21, start => 25585733, end => 25585733, allele_string => 'N/T', strand => -1});
is($p->validate_vf($vf), 1, 'validate_vf - lookup_ref -ve pass');
is($vf->allele_string, 'G/T', 'validate_vf - lookup_ref -ve seq OK');

$vf = get_vf({chr => 21, start => 25585733, end => 25585733, allele_string => 'N/-'});
is($p->validate_vf($vf), 1, 'validate_vf - lookup_ref deletion pass');
is($vf->allele_string, 'C/-', 'validate_vf - lookup_ref deletion seq OK');

$vf = get_vf({allele_string => '-/T', start => 2, end => 1});
is($p->validate_vf($vf), 1, 'validate_vf - lookup_ref insertion');
is($vf->allele_string, '-/T', 'validate_vf - lookup_ref insertion seq');

$p->{lookup_ref} = 0;


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

$p->file($test_cfg->{windows_vcf});
is($p->detect_format, 'vcf', 'detect_format - VCF with Windows newline');

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

$p->file($test_cfg->create_input_file('21   25587759 25587759     C/A  +  test'));
is($p->detect_format, 'ensembl', 'detect_format - multiple spaces delimiter');

$p->file($test_cfg->create_input_file('21 25587759 25587759 C/A/G + test'));
is($p->detect_format, 'ensembl', 'detect_format - VEP_input multiple');

$p->file($test_cfg->create_input_file('21 25587759 25587769 DUP + test'));
is($p->detect_format, 'ensembl', 'detect_format - VEP_input SV');

$p->file($test_cfg->create_input_file('21:25587759-25587759:1/A'));
is($p->detect_format, 'region', 'detect_format - region');

$p->file($test_cfg->create_input_file('21:25587759-25587759:1/DUP'));
is($p->detect_format, 'region', 'detect_format - region SV');

$p->file($test_cfg->create_input_file('chr21:25587759-25587759:1/A'));
is($p->detect_format, 'region', 'detect_format - region chr');

$p->file($test_cfg->create_input_file('CHR_21:25587759-25587759:1/A'));
is($p->detect_format, 'region', 'detect_format - region CHR_');

$p->file($test_cfg->create_input_file('21:25587759-25587759/A'));
is($p->detect_format, 'region', 'detect_format - region no strand');

$p->file($test_cfg->create_input_file('21:25587759-25587759:+1/A'));
is($p->detect_format, 'region', 'detect_format - region +ve strand');

$p->file($test_cfg->create_input_file('21:25587759-25587759:-1/A'));
is($p->detect_format, 'region', 'detect_format - region -ve strand');

$p->file($test_cfg->create_input_file('21:25587759-25587759:1:A'));
is($p->detect_format, 'region', 'detect_format - region : separator');

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

# detect format from STDIN dies
$tmp = '21 25587759 25587759 C/A + test';
open STDIN, '<', \$tmp;
$p->file('STDIN');
throws_ok {$p->detect_format} qr/Cannot detect format from STDIN/, 'fail detect format from STDIN';

# new with format configured
$p = Bio::EnsEMBL::VEP::Parser->new({config => $cfg, file => $test_cfg->{test_vcf}, format => 'vcf'});
is(ref($p), 'Bio::EnsEMBL::VEP::Parser::VCF', 'new with explicit format');

# new with invalid format
throws_ok {
  Bio::EnsEMBL::VEP::Parser->new({config => $cfg, file => $test_cfg->{test_vcf}, format => 'foo'})
} qr/Unknown or unsupported input format/, 'new with unknown format';

# new with format detection
$p = Bio::EnsEMBL::VEP::Parser->new({config => $cfg, file => $test_cfg->{test_vcf}, format => 'guess'});
is(ref($p), 'Bio::EnsEMBL::VEP::Parser::VCF', 'new with format detection');



## DATABASE TESTS
#################

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 11 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens') if $can_use_db;
  
  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$cfg_hash,
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    warning_file => 'STDERR',
  });

  $p = Bio::EnsEMBL::VEP::Parser->new({config => $cfg, file => $test_cfg->{test_vcf}});

  # contig to toplevel
  my $vf =  get_vf({allele_string => '-/C', chr => 'AP000235.3', start => 100, end => 99});
  is($p->validate_vf($vf), 1, 'DB - validate_vf - successful transform insertion 1');
  is($vf->{chr}, 21, 'DB - validate_vf - successful transform insertion 2');
  is($vf->{start}, 25043768, 'DB - validate_vf - successful transform insertion 3');

  $vf = get_vf({allele_string => 'A/C', chr => 'AP000235.3'});
  is($p->validate_vf($vf), 1, 'DB - validate_vf - successful transform 1');
  is($vf->{chr}, 21, 'DB - validate_vf - successful transform 2');
  is($vf->{start}, 25043669, 'DB - validate_vf - successful transform 3');

  # warning_msg prints to STDERR
  no warnings 'once';
  open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

  close STDERR;
  my $tmp;
  open STDERR, '>', \$tmp;

  is($p->validate_vf(get_vf({allele_string => 'A/C', chr => 2})), 0, 'DB - validate_vf - chromosome not in valid list 1');
  ok($tmp =~ /Could not fetch slice for chromosome/, 'DB - validate_vf - chromosome not in valid list 2');

  is($p->validate_vf(get_vf({allele_string => 'A/C', chr => 21, start => 125000001, end => 125000001})), 0, 'DB - validate_vf - unsuccessful transform to toplevel 1');
  ok($tmp =~ /could not transform to toplevel/, 'DB - validate_vf - unsuccessful transform to toplevel 2');

  # restore STDERR
  open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";


  # map to/from LRG
  $vf = get_vf({chr => '21', start => 43774213, end =>43774213, allele_string => 'T/C'});
  my $mapped = $p->map_to_lrg([$vf]);
  is(scalar @$mapped, 2, 'map_to_lrg - count');
  is($mapped->[1]->{chr}, 'LRG_485', 'map_to_lrg - check chr');

  $vf = get_vf({chr => 'LRG_485', start => 7166, end =>7166, allele_string => 'T/C'});
  $mapped = $p->map_to_lrg([$vf]);
  is(scalar @$mapped, 2, 'map_to_lrg - from LRG - count');
  is($mapped->[1]->{chr}, '21', 'map_to_lrg - from LRG - check chr');
}


done_testing();


sub get_vf {
  my $hashref = shift;

  $hashref->{$_} ||= 1 for qw(chr start end strand);

  return Bio::EnsEMBL::Variation::VariationFeature->new_fast($hashref);
}
