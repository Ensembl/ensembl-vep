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

use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::StructuralVariationFeature;

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

my $cfg_hash = $test_cfg->base_testing_cfg;

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::OutputFactory::VCF');

use_ok('Bio::EnsEMBL::VEP::Config');
use_ok('Bio::EnsEMBL::VEP::Runner');
my $cfg = Bio::EnsEMBL::VEP::Config->new($cfg_hash);

my $of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $cfg, header_info => $test_cfg->{header_info}});

is(ref($of), 'Bio::EnsEMBL::VEP::OutputFactory::VCF', 'check class');



## METHOD TESTS
###############

is_deeply(
  $of->headers,
  [
    '##fileformat=VCFv4.1',
    '##VEP="v1" time="test"',
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|custom_test">',
    '##INFO=<ID=custom_test,Number=.,Type=String,Description="test.vcf.gz (overlap)">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ],
  'headers'
);

my $headers = get_runner({plugin => ['TestPlugin'], quiet => 1, input_file => $test_cfg->{test_vcf}, vcf => 1})->get_OutputFactory->headers;
is_deeply(
  [$headers->[-3], $headers->[-2], $headers->[-1]],
  [
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|test">',
    '##test=header',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096"
  ],
  'headers - plugin'
);

is_deeply(
  $of->fields,
  [
    'Allele',
    'Consequence',
    'IMPACT',
    'SYMBOL',
    'Gene',
    'Feature_type',
    'Feature',
    'BIOTYPE',
    'EXON',
    'INTRON',
    'HGVSc',
    'HGVSp',
    'cDNA_position',
    'CDS_position',
    'Protein_position',
    'Amino_acids',
    'Codons',
    'Existing_variation',
    'DISTANCE',
    'STRAND',
    'FLAGS',
    'custom_test'
  ],
  'fields - default'
);

$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({sift => 'b'})
});
is_deeply(
  $of->fields,
  [
    'Allele',
    'Consequence',
    'IMPACT',
    'SYMBOL',
    'Gene',
    'Feature_type',
    'Feature',
    'BIOTYPE',
    'EXON',
    'INTRON',
    'HGVSc',
    'HGVSp',
    'cDNA_position',
    'CDS_position',
    'Protein_position',
    'Amino_acids',
    'Codons',
    'Existing_variation',
    'DISTANCE',
    'STRAND',
    'FLAGS',
    'SIFT',
  ],
  'fields - flag on'
);

$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({
  config => Bio::EnsEMBL::VEP::Config->new({fields => 'Allele,Consequence'})
});
is_deeply(
  $of->fields,
  [
    'Allele',
    'Consequence',
  ],
  'fields - user-defined'
);

$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $cfg});



## output_hash_to_vcf_info_chunk
################################

is(
  $of->output_hash_to_vcf_info_chunk({}),
  ('|' x 20),
  'output_hash_to_vcf_info_chunk - empty'
);

is(
  $of->output_hash_to_vcf_info_chunk({Allele => 'A'}),
  'A'.('|' x 20),
  'output_hash_to_vcf_info_chunk'
);

is(
  $of->output_hash_to_vcf_info_chunk({Allele => 'A'}, -1),
  'T'.('|' x 20),
  'output_hash_to_vcf_info_chunk - allele revcomped'
);

is(
  $of->output_hash_to_vcf_info_chunk({Consequence => '-'}),
  ('|' x 20),
  'output_hash_to_vcf_info_chunk - "-" erased'
);

is(
  $of->output_hash_to_vcf_info_chunk({Allele => '-'}),
  '-'.('|' x 20),
  'output_hash_to_vcf_info_chunk - "-" intact for Allele'
);

is(
  $of->output_hash_to_vcf_info_chunk({Allele => ';'}),
  '%3B'.('|' x 20),
  'output_hash_to_vcf_info_chunk - ";" uri encoded'
);

## VariationFeature_to_VCF_record
#################################

is_deeply(
  $of->VariationFeature_to_VCF_record(get_vf({allele_string => 'A/G'})),
  [1, 1, '.', 'A', 'G', '.', '.', '.'],
  'VariationFeature_to_VCF_record'
);

is_deeply(
  $of->VariationFeature_to_VCF_record(get_vf({allele_string => 'A/G', strand => -1})),
  [1, 1, '.', 'T', 'C', '.', '.', '.'],
  'VariationFeature_to_VCF_record - rev strand'
);

is_deeply(
  $of->VariationFeature_to_VCF_record(get_vf({allele_string => 'A/G/T'})),
  [1, 1, '.', 'A', 'G,T', '.', '.', '.'],
  'VariationFeature_to_VCF_record - multiple alts'
);

is_deeply(
  $of->VariationFeature_to_VCF_record(get_vf({allele_string => 'AG/CT'})),
  [1, 1, '.', 'AG', 'CT', '.', '.', '.'],
  'VariationFeature_to_VCF_record - balanced non-SNP'
);

is_deeply(
  $of->VariationFeature_to_VCF_record(get_vf({allele_string => 'A/-', start => 2, end => 2})),
  [1, 1, '.', 'NA', 'N', '.', '.', '.'],
  'VariationFeature_to_VCF_record - deletion - no seq lookup'
);

is_deeply(
  $of->VariationFeature_to_VCF_record(get_vf({allele_string => '-/A', start => 2})),
  [1, 1, '.', 'N', 'NA', '.', '.', '.'],
  'VariationFeature_to_VCF_record - insertion - no seq lookup'
);

is_deeply(
  $of->VariationFeature_to_VCF_record(get_vf({allele_string => 'A/-/G', start => 2, end => 2})),
  [1, 1, '.', 'NA', 'N,NG', '.', '.', '.'],
  'VariationFeature_to_VCF_record - mixed - no seq lookup'
);

is_deeply(
  $of->VariationFeature_to_VCF_record(
    Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
      class_SO_term => 'deletion',
      start => 11,
      end => 20,
      chr => 1,
    })
  ),
  [1, 10, '.', 'N', '<DEL>', '.', '.', 'END=20'],
  'VariationFeature_to_VCF_record - SV'
);

# test with sequence fetch
is_deeply(
  get_runner({
    input_file => $test_cfg->{test_vcf},
    dir => $test_cfg->{cache_root_dir},
    vcf => 1,
  })->get_OutputFactory->VariationFeature_to_VCF_record(
    get_vf({
      allele_string => 'A/-',
      start => 25585733,
      end => 25585733,
      chr => 21,
    })
  ),
  [21, 25585732, '.', 'GA', 'G', '.', '.', '.'],
  'VariationFeature_to_VCF_record - deletion - seq lookup'
);


## get_all_lines_by_InputBuffer
###############################

my $runner = get_runner({input_file => $test_cfg->{test_vcf}});
my $ib = $runner->get_InputBuffer;
$of = $runner->get_OutputFactory;

my @lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(scalar @lines, scalar @{$ib->buffer}, 'get_all_lines_by_InputBuffer - count');

is(
  $lines[0],
  "21\t25585733\trs142513484\tC\tT\t.\t.\t".
  'CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|'.
  "\tGT\t0|0",
  'get_all_lines_by_InputBuffer - check first'
);

is(
  $lines[-1],
  "21\t25982445\trs141331202\tC\tT\t.\t.\t".
  'CSQ=T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000346798||||||1157|1123|375|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000348990||||||1045|898|300|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000354192||||||857|730|244|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000357903||||||1233|1066|356|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000358918||||||1170|1123|375|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000359726||||||985|793|265|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000415997||||||326|328|110|V/I|Gtt/Att|||-1|cds_start_NF&cds_end_NF,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000439274||||||989|955|319|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000440126||||||1317|1051|351|V/I|Gtt/Att|||-1|,'.
  'T|missense_variant|MODERATE||ENSG00000142192|Transcript|ENST00000448850||||||830|832|278|V/I|Gtt/Att|||-1|cds_start_NF'.
  "\tGT\t0|0",
  'get_all_lines_by_InputBuffer - check last'
);

$ib = get_runner({
  input_file => $test_cfg->{test_vcf},
  everything => 1,
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;

$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});

@lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(
  $lines[0],
  "21\t25585733\trs142513484\tC\tT\t.\t.\t".
  'CSQ=T|3_prime_UTR_variant|MODIFIER|MRPL39|ENSG00000154719|Transcript|ENST00000307301|protein_coding|'.
  '11/11||ENST00000307301.11:c.*18G>A||1122|||||rs142513484||-1||SNV|HGNC|HGNC:14027|YES|5|A2|CCDS33522.1|'.
  'ENSP00000305682|Q9NYK5||UPI00001AEAC0||||||0.0010|0.0030|0.0014|0.0000|0.0000|0.0000|0.005|'.
  '0|4.119e-04|0.0004133|0.004681|0.000173|0|0|0|0|0||||||||,'.
  'T|missense_variant|MODERATE|MRPL39|ENSG00000154719|Transcript|ENST00000352957|protein_coding|'.
  '10/10||ENST00000352957.8:c.991G>A|ENSP00000284967.6:p.Ala331Thr|1033|991|331|A/T|Gca/Aca|rs142513484|'.
  '|-1||SNV|HGNC|HGNC:14027||1|P3|CCDS13573.1|ENSP00000284967|Q9NYK5||UPI00001AEE66||'.
  'tolerated_low_confidence(0.17)|benign(0.021)|||0.0010|0.0030|0.0014|0.0000|0.0000|0.0000|0.005|'.
  '0|4.119e-04|0.0004133|0.004681|0.000173|0|0|0|0|0||||||||,'.
  'T|upstream_gene_variant|MODIFIER|AP000223.42|ENSG00000260583|Transcript|ENST00000567517|antisense|'.
  '|||||||||rs142513484|2407|-1||SNV|Clone_based_vega_gene||YES|||||||||||||0.0010|0.0030|0.0014|0.0000|'.
  '0.0000|0.0000|0.005|0|4.119e-04|0.0004133|0.004681|0.000173|0|0|0|0|0||||||||,'.
  'T|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001963192|TF_binding_site||||||||||rs142513484|'.
  '|||SNV||||||||||||||||0.0010|0.0030|0.0014|0.0000|0.0000|0.0000|0.005|0|4.119e-04|0.0004133|'.
  '0.004681|0.000173|0|0|0|0|0||||||||'.
  "\tGT\t0|0",
  'get_all_lines_by_InputBuffer - everything'
);

# custom
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 1 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

  $runner = get_runner({
    input_file => $test_cfg->{test_vcf},
    custom => [$test_cfg->{custom_vcf}.',test,vcf,exact,,FOO'],
    output_format => 'vcf',
  });
  $of = $runner->get_OutputFactory;

  is(
    $of->get_all_lines_by_InputBuffer($runner->get_InputBuffer)->[0],
    "21\t25585733\trs142513484\tC\tT\t.\t.\t".
    "CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|||test1|BAR,".
    "T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|||test1|BAR,".
    "T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|||test1|BAR\tGT\t0|0",
    'get_all_lines_by_InputBuffer - custom'
  );
}

# test converting to VCF from different input
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 25585733 C/T 1)]),
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});

is_deeply(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585733\t21_25585733_C/T\tC\tT\t.\t.\t".
  'CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  "non-VCF input"
);

# test keep vs trash existing CSQ
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 . C T . . CSQ=foo;BAR=blah)]),
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});

is_deeply(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585733\t.\tC\tT\t.\t.\t".
  'BAR=blah;CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  "trash existing CSQ"
);

$of->{keep_csq} = 1;
is_deeply(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585733\t.\tC\tT\t.\t.\t".
  'CSQ=foo;BAR=blah;CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  "keep existing CSQ"
);

# different VCF info field
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 . C T . . .)]),
  dir => $test_cfg->{cache_root_dir},
  vcf_info_field => 'EFF',
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});

@lines = @{$of->get_all_lines_by_InputBuffer($ib)};
ok($lines[0] =~ /EFF/ && $lines[0] !~ /CSQ/, 'vcf_info_field');




## test getting stuff from input
################################

$runner = get_runner({
  input_file => $test_cfg->{test_vcf},
  dir => $test_cfg->{cache_root_dir},
  vcf => 1,
});

$of = $runner->get_OutputFactory;

is(
  $of->headers->[0],
  '##fileformat=VCFv4.1',
  'headers - from input 1'
);


is(
  $of->headers->[2],
  '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">',
  'headers - from input 2'
);

is(
  $of->headers->[3],
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096",
  'headers - from input 3'
);



## web_output
#############

my $tmp_file = $test_cfg->create_input_file();

$ib = get_runner({
  input_file => $test_cfg->{test_vcf},
  dir => $test_cfg->{cache_root_dir},
  web_output => $tmp_file
})->get_InputBuffer;

Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config})->get_all_lines_by_InputBuffer($ib);

open IN, $tmp_file;
@lines = <IN>;
close IN;

is(scalar @lines, 132, 'web_output - count lines');
is($lines[-1], "21\t25982445\t25982445\tC/T\t1\trs141331202\tmissense_variant\n", 'web_output - check last');

# test long allele string truncation
$ib = get_runner({
  input_file => $test_cfg->{test_vcf},
  dir => $test_cfg->{cache_root_dir},
  web_output => $tmp_file
})->get_InputBuffer;

$ib->buffer->[0]->{allele_string} = 'C/'.('T' x 100);

Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config})->get_all_lines_by_InputBuffer($ib);

open IN, $tmp_file;
@lines = <IN>;
close IN;

is($lines[0], "21\t25585733\t25585733\tC/100BP_SEQ\t1\trs142513484\tmissense_variant\n", 'web_output - allele truncation');



# done
done_testing();


sub get_vf {
  my $hashref = shift;

  $hashref->{$_} ||= 1 for qw(chr start end strand);

  return Bio::EnsEMBL::Variation::VariationFeature->new_fast($hashref);
}

sub get_runner {
  my $tmp_cfg = shift;

  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir},
    %$tmp_cfg,
    output_format => 'vcf'
  });

  $runner->init;

  my $ib = $runner->get_InputBuffer;
  $ib->next();
  $_->annotate_InputBuffer($ib) for @{$runner->get_all_AnnotationSources};
  $ib->finish_annotation();

  return $runner;
}
