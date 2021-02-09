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
use_ok('Bio::EnsEMBL::VEP::Parser::VCF');

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($base_testing_cfg);
ok($cfg, 'get new config object');

my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
ok($p, 'new is defined');

is(ref($p), 'Bio::EnsEMBL::VEP::Parser::VCF', 'check class');



## METHOD TESTS
###############

is(ref($p->parser), 'Bio::EnsEMBL::IO::Parser::VCF4', 'parser');

is_deeply(
  $p->headers,
  [
    '##fileformat=VCFv4.1',
    '##contig=<ID=21,assembly=GCF_000001405.26,length=46709983>',
    '##contig=<ID=22,assembly=GCF_000001405.26,length=50818468>',
    '##ALT=<ID=CNV,Description="Copy Number Polymorphism">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '#'.join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT HG00096)),
  ],
  'headers'
);

$vf = $p->next;
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'rs142513484',
  'map_weight' => 1,
  'allele_string' => 'C/T',
  'end' => 25585733,
  'start' => 25585733,
  'seq_region_end' => 25585733,
  'seq_region_start' => 25585733
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'next');

is(ref($p->next), 'Bio::EnsEMBL::Variation::VariationFeature', 'next again');

# check not shorting out on non-variant lines
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([
    [qw(21 25587759 test1 C A . . .)],
    [qw(21 25587760 test2 C . . . .)],
    [qw(21 25587761 test3 C . . . .)],
    [qw(21 25587762 test4 C A . . .)]
  ]),
  valid_chromosomes => [21]
});

is($p->next->variation_name, 'test1', 'no shorting out 1');
is($p->next->variation_name, 'test4', 'no shorting out 2');


## FORMAT TESTS
###############

# disable warnings for using "," in a qw()
no warnings 'qw';

# deletion
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test AC A . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/-',
  'end' => 25587760,
  'start' => 25587760,
  'seq_region_end' => 25587760,
  'seq_region_start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'deletion');

# insertion
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test A AC . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => '-/C',
  'end' => 25587759,
  'start' => 25587760,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'insertion');

# deletion with SVTYPE field
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test AC A . . SVTYPE=DEL)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/-',
  'end' => 25587760,
  'start' => 25587760,
  'seq_region_end' => 25587760,
  'seq_region_start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'deletion ignore SVTYPE');

# deletion with SVTYPE field
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test A AC . . SVLEN=1)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => '-/C',
  'end' => 25587759,
  'start' => 25587760,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'insertion ignore SVLEN');

# multiple alts
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test A C,G . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/C/G',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'multiple alts');

# mixed types - different first base
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test A C,GG . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/C/GG',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'mixed types - different first base');

# mixed types - different first base
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G GC,GT . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => '-/C/T',
  'end' => 25587759,
  'start' => 25587760,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'mixed types - same first base');

# stubby - nothing after ALT
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G C)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'G/C',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'stubby');

# non-variant
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, allow_non_variant => 1}),
  file => $test_cfg->create_input_file([qw(21 25587759 test G . . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'G',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759,
  'non_variant' => 1,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'non-variant');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, allow_non_variant => 0}),
  file => $test_cfg->create_input_file([qw(21 25587759 test G . . . .)]),
  valid_chromosomes => [21]
})->next();
is($vf, undef, 'non-variant without allow_non_variant');

# *-type as produced by GATK
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G C,* . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'G/C/*',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type 1');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G C,<DEL:*> . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'G/C/<DEL:*>',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type 2 ');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test GC G,* . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});

is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/-/*',
  'end' => 25587760,
  'start' => 25587760,
  'seq_region_end' => 25587760,
  'seq_region_start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type with deletion');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G GC,* . . .)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});

is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => '-/C/*',
  'end' => 25587759,
  'start' => 25587760,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type with insertion');

# minimal
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, minimal => 1}),
  file => $test_cfg->create_input_file([qw(21 25587758 test CAT CCT . . .)]),
  valid_chromosomes => [21],
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'minimised' => 1,
  'original_allele_string' => 'CAT/CCT',
  'original_end' => 25587760,
  'end' => 25587759,
  'seq_region_end' => 25587759,
  'original_start' => '25587758',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/C',
  'start' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'minimal');


$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, minimal => 1}),
  file => $test_cfg->create_input_file([qw(21 25587758 test C T,CAA . . .)]),
  valid_chromosomes => [21],
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/T/CAA',
  'end' => 25587758,
  'start' => 25587758,
  'seq_region_end' => 25587758,
  'seq_region_start' => 25587758
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'minimal - multiple alts pass through here');


# basic SV coord tests
$expected = bless( {
  'outer_end' => 25587769,
  'chr' => '21',
  'inner_end' => 25587769,
  'outer_start' => 25587759,
  'end' => 25587769,
  'seq_region_end' => 25587769,
  'inner_start' => 25587759,
  'strand' => 1,
  'class_SO_term' => 'duplication',
  'variation_name' => 'sv_dup',
  'start' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' );

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . SVTYPE=DUP;END=25587769)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected, 'StructuralVariationFeature 1');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T . . . SVTYPE=DUP;END=25587769)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected, 'StructuralVariationFeature 2');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . END=25587769)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected , 'StructuralVariationFeature no SVTYPE');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . SVLEN=11)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected , 'StructuralVariationFeature SVLEN');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . SVLEN=11;CIPOS=-3,2;CIEND=-4,5)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'outer_end' => 25587774,
  'chr' => '21',
  'inner_end' => 25587765,
  'outer_start' => 25587756,
  'end' => 25587769,
  'seq_region_end' => 25587769,
  'inner_start' => 25587761,
  'strand' => 1,
  'class_SO_term' => 'duplication',
  'variation_name' => 'sv_dup',
  'start' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' ) , 'StructuralVariationFeature fuzzy');

no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DEL> . . .)]),
  valid_chromosomes => [21]
})->next();
ok($tmp =~ /VCF line.+looks incomplete/, 'StructuralVariationFeature del without end or length');

open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

# Test if incomplete variant (DEL) is skipped
open(SAVE, ">&STDERR") or die "Can't save STDERR\n";
close STDERR;
open STDERR, '>', \$tmp;
my $vf_del = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([
    [qw(21 25587758 sv_del T <DEL> . . .)],
    [qw(21 25587759 test A C . . .)]
    ]),
  valid_chromosomes => [21]
});

my $sv = $vf_del->next();
ok($tmp =~ /VCF line.+looks incomplete/, 'StructuralVariationFeature del without end or length (2 variants)');
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

delete($sv->{adaptor});
delete($sv->{_line});

is_deeply($sv, bless( {
  'outer_end' => 25587759,
  'chr' => '21',
  'inner_end' => 25587759,
  'outer_start' => 25587759,
  'end' => 25587759,
  'vep_skip' => 1,
  'seq_region_end' => 25587759,
  'inner_start' => 25587759,
  'strand' => 1,
  'class_SO_term' => 'deletion',
  'variation_name' => 'sv_del',
  'start' => 25587759,
  'seq_region_start' => 25587759,
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' ) , 'StructuralVariationFeature - skipping incomplete variant');

my $snv = $vf_del->next();
delete($snv->{adaptor});
delete($snv->{_line});

is_deeply($snv, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/C',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'VariationFeature - variant not skipped');

## test max SV length
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n";
close STDERR;
open STDERR, '>', \$tmp;

my $lvf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1, max_sv_size => 1000, warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . SVLEN=10001;CIPOS=-3,2;CIEND=-4,5)]),
  valid_chromosomes => [21]
})->next();
delete($lvf->{adaptor}); delete($lvf->{_line});

is_deeply($lvf, bless( {
  'outer_end' => 25597764,
  'chr' => '21',
  'inner_end' => 25597755,
  'outer_start' => 25587756,
  'end' => 25597759,
  'vep_skip' => 1,
  'seq_region_end' => 25597759,
  'inner_start' => 25587761,
  'strand' => 1,
  'class_SO_term' => 'duplication',
  'variation_name' => 'sv_dup',
  'start' => 25587759,
  'seq_region_start' => 25587759,
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' ) , 'StructuralVariationFeature - longer than specified maximum');

open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

## test complex SV
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n";
close STDERR;
open STDERR, '>', \$tmp;

my $cvf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1, max_sv_size => 1000, warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(1 774569 gnomAD_v2_CPX_1_1 N	<CPX> 1 PASS END=828435;SVTYPE=CPX;CHR2=1;SVLEN=53959)]),
  valid_chromosomes => [1]
})->next();
delete($cvf->{adaptor}); delete($cvf->{_line});

is_deeply($cvf, bless( {
                 'outer_end' => '828435',
                 'chr' => '1',
                 'inner_end' => '828435',
                 'outer_start' => '774570',
                 'end' => 828435,
                 'vep_skip' => 1,
                 'seq_region_start' => 774570,
                 'inner_start' => '774570',
                 'strand' => 1,
                 'seq_region_end' => 828435,
                 'class_SO_term' => 'CPX',
                 'variation_name' => 'gnomAD_v2_CPX_1_1',
                 'start' => 774570
               },
                'Bio::EnsEMBL::Variation::StructuralVariationFeature' ) , 'StructuralVariationFeature - CPX skipped');


ok($tmp =~ /variant CPX is of a non-supported type/, 'StructuralVariationFeature - skip CPX warning');

open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";



my $cnv_vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1,  warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(1 774569 gnomAD_v2_DEL_1_1 N <CN0> 1 PASS END=828435;SVTYPE=DUP;CHR2=1;SVLEN=53959)]),
  valid_chromosomes => [1]
})->next();
delete($cnv_vf->{adaptor}); delete($cnv_vf->{_line});

is_deeply($cnv_vf, bless( {
                 'outer_end' => '828435',
                 'chr' => '1',
                 'inner_end' => '828435',
                 'outer_start' => '774570',
                 'end' => 828435,
                 'seq_region_start' => 774570,
                 'inner_start' => '774570',
                 'strand' => 1,
                 'seq_region_end' => 828435,
                 'class_SO_term' => 'duplication',
                 'variation_name' => 'gnomAD_v2_DEL_1_1',
                 'start' => 774570
               },
                'Bio::EnsEMBL::Variation::StructuralVariationFeature' ) , 'StructuralVariationFeature - CNV with only duplication allele');

$cnv_vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1,  warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(1 774569 gnomAD_v2_DEL_1_1 N <CN0>,<CN2> 1 PASS END=828435;SVTYPE=CNV;CHR2=1;SVLEN=53959)]),
  valid_chromosomes => [1]
})->next();
delete($cnv_vf->{adaptor}); delete($cnv_vf->{_line});

is_deeply($cnv_vf, bless( {
                 'outer_end' => '828435',
                 'chr' => '1',
                 'inner_end' => '828435',
                 'outer_start' => '774570',
                 'end' => 828435,
                 'seq_region_start' => 774570,
                 'inner_start' => '774570',
                 'strand' => 1,
                 'seq_region_end' => 828435,
                 'class_SO_term' => 'copy_number_variation',
                 'variation_name' => 'gnomAD_v2_DEL_1_1',
                 'start' => 774570
               },
                'Bio::EnsEMBL::Variation::StructuralVariationFeature' ) , 'StructuralVariationFeature - CNV with multiple alleles');


## Currently regard fusions across different chromosomes as single breakpoint
my $bnd_vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1,  warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(2	68914092	BND00001121	A	]1:37938377]A	.	PASS	SVTYPE=BND;CHR2=1;END=37938377   )]),
  valid_chromosomes => [1,2]
})->next();
ok($bnd_vf->end() == 68914093, 'StructuralVariationFeature - cross-chromosome BND');

my $bnd2_vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1,  warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(2    68914092        BND00001121     A       ]2:68920000]A	.	PASS   SVTYPE=BND;CHR2=2;END=68920000   )]),
  valid_chromosomes => [1,2]
})->next();
ok($bnd2_vf->end() == 68920000, 'StructuralVariationFeature - within-chromosome BND');




## OTHER TESTS
##############

# GP flag
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1}),
  file => $test_cfg->create_input_file([qw(21 25587759 test A G . . GP=21:25586000)]),
  valid_chromosomes => [21]
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => 25586000,
  'start' => 25586000,
  'seq_region_end' => 25586000,
  'seq_region_start' => 25586000
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'gp');

# GP flag not found
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, gp => 1, warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(21 25587759 test A G . . .)]),
  valid_chromosomes => [21]
})->next();
ok($tmp =~ /No GP flag found in INFO column/, 'gp - not found');

open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";


# individual data
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, allow_non_variant => 1, individual => 'all'}),
  file => $test_cfg->create_input_file([
    ['##fileformat=VCFv4.1'],
    [qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT dave barry jeff)],
    [qw(21 25587759 indtest A G . . . GT 0|1 1/1 0/0)],
  ]),
  valid_chromosomes => [21]
});

$vf = $p->next;
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'indtest',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759,
  'genotype' => ['A', 'G'],
  'individual' => 'dave',
  'phased' => 1,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'individual 1');

$vf = $p->next;
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'indtest',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759,
  'genotype' => ['G', 'G'],
  'individual' => 'barry',
  'phased' => 0,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'individual 2');

$vf = $p->next;
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'indtest',
  'map_weight' => 1,
  'allele_string' => 'A',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759,
  'genotype' => ['A', 'A'],
  'individual' => 'jeff',
  'phased' => 0,
  'non_variant' => 1,
  'hom_ref' => 1,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'individual 3');


$p = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, allow_non_variant => 1, process_ref_homs => 1, individual => 'jeff'}),
  file => $test_cfg->create_input_file([
    ['##fileformat=VCFv4.1'],
    [qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT jeff)],
    [qw(21 25587759 indtest A G . . . GT 0/0)],
  ]),
  valid_chromosomes => [21]
});

$vf = $p->next;
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'indtest',
  'map_weight' => 1,
  'allele_string' => 'A/A',
  'end' => 25587759,
  'start' => 25587759,
  'seq_region_end' => 25587759,
  'seq_region_start' => 25587759,
  'genotype' => ['A', 'A'],
  'individual' => 'jeff',
  'phased' => 0,
  'hom_ref' => 1,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'individual - process_ref_homs');


# test situation where a line gets passed over for having no valid ALT genotypes
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, individual => 'jeff'}),
  file => $test_cfg->create_input_file([
    ['##fileformat=VCFv4.1'],
    [qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT jeff)],
    [qw(21 25587759 indtest1 A G . . . GT ./.)],
    [qw(21 25587760 indtest2 C T . . . GT 1/1)],
  ]),
  valid_chromosomes => [21]
});

$vf = $p->next;
is($vf->{variation_name}, 'indtest2', 'individual - pass over entries with no valid GTs');


# test alleles being mapped properly to ensembl types
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$base_testing_cfg, individual => 'jeff'}),
  file => $test_cfg->create_input_file([
    ['##fileformat=VCFv4.1'],
    [qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT jeff)],
    [qw(21 25587759 indtest1 A AT . . . GT 1/1)],
    [qw(21 25587760 indtest2 GC G . . . GT 1/1)],
  ]),
  valid_chromosomes => [21]
});

$vf = $p->next;
is($vf->{allele_string}, '-/T', 'individual - insertion alleles mapped');

$vf = $p->next;
is($vf->{allele_string}, 'C/-', 'individual - deletion alleles mapped');





# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$cfg->param('warning_file', 'STDERR');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([
    [qw(21 foo . C A . . .)],
    [qw(21 25587759 . C A . . .)],
  ]),
  valid_chromosomes => [21]
})->next();

is($vf->{start}, 25587759, 'skip VF that fails validation');

$cfg->param('dont_skip', 1);

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([
    [qw(21 foo . C A . . .)],
    [qw(21 25587759 . C A . . .)],
  ]),
  valid_chromosomes => [21]
})->next();

is($vf->{start}, 'foo', 'dont skip VF that fails validation with dont_skip');

# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";


# done
done_testing();
