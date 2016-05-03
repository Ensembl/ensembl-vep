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

my ($vf, $tmp, $expected);

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

is_deeply(
  $p->headers,
  [
    ['fileformat', 'VCFv4.1'],
    [
      'header', [
        'CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO',
        'FORMAT',
        'HG00096'
      ],
    ],
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
  'start' => 25585733
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'next');

is(ref($p->next), 'Bio::EnsEMBL::Variation::VariationFeature', 'next again');



## FORMAT TESTS
###############

# disable warnings for using "," in a qw()
no warnings 'qw';

# deletion
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test AC A . . .)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'C/-',
  'end' => 25587760,
  'start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'deletion');

# insertion
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test A AC . . .)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => '-/C',
  'end' => 25587759,
  'start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'insertion');

# multiple alts
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test A C,G . . .)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/C/G',
  'end' => 25587759,
  'start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'multiple alts');

# mixed types - different first base
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test A C,GG . . .)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/C/GG',
  'end' => 25587759,
  'start' => 25587759
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'mixed types - different first base');

# mixed types - different first base
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G GC,GT . . .)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => '-/C/T',
  'end' => 25587759,
  'start' => 25587760
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'mixed types - same first base');

# non-variant
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({allow_non_variant => 1}),
  file => $test_cfg->create_input_file([qw(21 25587759 test G . . . .)])
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
  'non_variant' => 1,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'non-variant');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({allow_non_variant => 0}),
  file => $test_cfg->create_input_file([qw(21 25587759 test G . . . .)])
})->next();
is($vf, undef, 'non-variant without allow_non_variant');

# *-type as produced by GATK
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G C,* . . .)])
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
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type 1');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G C,<DEL:*> . . .)])
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
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type 2 ');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test GC G,* . . .)])
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
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type with deletion');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587759 test G GC,* . . .)])
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
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), '*-type with insertion');


# basic SV coord tests
$expected = bless( {
  'outer_end' => 25587769,
  'chr' => '21',
  'inner_end' => 25587769,
  'outer_start' => 25587759,
  'end' => 25587769,
  'inner_start' => 25587759,
  'strand' => 1,
  'class_SO_term' => 'duplication',
  'variation_name' => 'sv_dup',
  'start' => 25587759
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' );

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . SVTYPE=DUP;END=25587769)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected, 'StructuralVariationFeature 1');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T . . . SVTYPE=DUP;END=25587769)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected, 'StructuralVariationFeature 2');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . END=25587769)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected , 'StructuralVariationFeature no SVTYPE');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . SVLEN=11)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, $expected , 'StructuralVariationFeature SVLEN');

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $cfg,
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DUP> . . SVLEN=11;CIPOS=-3,2;CIEND=-4,5)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'outer_end' => 25587774,
  'chr' => '21',
  'inner_end' => 25587765,
  'outer_start' => 25587756,
  'end' => 25587769,
  'inner_start' => 25587761,
  'strand' => 1,
  'class_SO_term' => 'duplication',
  'variation_name' => 'sv_dup',
  'start' => 25587759
}, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' ) , 'StructuralVariationFeature fuzzy');

no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({gp => 1, warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(21 25587758 sv_dup T <DEL> . . .)])
})->next();
ok($tmp =~ /VCF line.+looks incomplete/, 'StructuralVariationFeature del without end or length');

open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";



## OTHER TESTS
##############

# GP flag
$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({gp => 1}),
  file => $test_cfg->create_input_file([qw(21 25587759 test A G . . GP=21:25586000)])
})->next();
delete($vf->{adaptor}); delete($vf->{_line});
is_deeply($vf, bless( {
  'chr' => '21',
  'strand' => 1,
  'variation_name' => 'test',
  'map_weight' => 1,
  'allele_string' => 'A/G',
  'end' => 25586000,
  'start' => 25586000
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'gp');

# GP flag not found
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$vf = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({gp => 1, warning_file => 'STDERR'}),
  file => $test_cfg->create_input_file([qw(21 25587759 test A G . . .)])
})->next();
ok($tmp =~ /No GP flag found in INFO column/, 'gp - not found');

open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";


# individual data
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({allow_non_variant => 1, individual => 'all'}),
  file => $test_cfg->create_input_file([
    ['##fileformat=VCFv4.1'],
    [qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT dave barry jeff)],
    [qw(21 25587759 indtest A G . . . GT 0|1 1/1 0/0)],
  ])
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
  'genotype' => ['A', 'A'],
  'individual' => 'jeff',
  'phased' => 0,
  'non_variant' => 1,
  'hom_ref' => 1,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'individual 3');


$p = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({allow_non_variant => 1, process_ref_homs => 1, individual => 'jeff'}),
  file => $test_cfg->create_input_file([
    ['##fileformat=VCFv4.1'],
    [qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT jeff)],
    [qw(21 25587759 indtest A G . . . GT 0/0)],
  ])
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
  'genotype' => ['A', 'A'],
  'individual' => 'jeff',
  'phased' => 0,
  'hom_ref' => 1,
}, 'Bio::EnsEMBL::Variation::VariationFeature' ), 'individual - process_ref_homs');


# done
done_testing();