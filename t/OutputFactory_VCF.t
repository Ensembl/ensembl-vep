# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

$headers = get_runner({refseq => 1, fasta => $test_cfg->{fasta}, quiet => 1, input_file => $test_cfg->{test_vcf}, vcf => 1})->get_OutputFactory->headers;
is(
  $headers->[-2],
  '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|GIVEN_REF|USED_REF|BAM_EDIT">',
  'headers - BAM_EDIT'
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
  config => Bio::EnsEMBL::VEP::Config->new({merged => 1, custom => 1, database => 0})
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
    'REFSEQ_MATCH',
    'SOURCE',
  ],
  'fields - --merged and --custom dont duplicate SOURCE'
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

is(
  $of->output_hash_to_vcf_info_chunk({Allele => '|'}),
  '&'.('|' x 20),
  'output_hash_to_vcf_info_chunk - "|" converted'
);

is(
  $of->output_hash_to_vcf_info_chunk({Allele => ','}),
  '&'.('|' x 20),
  'output_hash_to_vcf_info_chunk - "," converted'
);

is(
  $of->output_hash_to_vcf_info_chunk({Allele => 'A  G'}),
  'A_G'.('|' x 20),
  'output_hash_to_vcf_info_chunk - whitespace converted'
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


# check incomplete VCF entry gets filled out properly
$ib = get_runner({
  input_data => qq{21\t25585733\trs142513484\tC\tT},
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;

$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});
@lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(
  $lines[0],
  "21\t25585733\trs142513484\tC\tT\t.\t.\t".
  'CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  'get_all_lines_by_InputBuffer - incomplete VCF entry filled out'
);

# check character conversion
$ib->buffer->[0]->get_all_TranscriptVariations->[0]->transcript->{stable_id} = 'ENST00,00  03|07;301';
@lines = @{$of->get_all_lines_by_InputBuffer($ib)};
is(
  $lines[0],
  "21\t25585733\trs142513484\tC\tT\t.\t.\t".
  'CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00&00_03&07%3B301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  'get_all_lines_by_InputBuffer - invalid character conversion'
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
  'ENSP00000305682|Q9NYK5||UPI00001AEAC0||||||0.0010|0.003|0.0014|0|0|0|0.004998|0|0.0003478|0.004643|'.
  '0.0003236|0|0|0|1.886e-05|0|0|0.004998|AA||||||||,'.
  'T|missense_variant|MODERATE|MRPL39|ENSG00000154719|Transcript|ENST00000352957|protein_coding|'.
  '10/10||ENST00000352957.8:c.991G>A|ENSP00000284967.6:p.Ala331Thr|1033|991|331|A/T|Gca/Aca|rs142513484|'.
  '|-1||SNV|HGNC|HGNC:14027||1|P3|CCDS13573.1|ENSP00000284967|Q9NYK5||UPI00001AEE66||tolerated_low_confidence(0.17)|'.
  'benign(0.021)|||0.0010|0.003|0.0014|0|0|0|0.004998|0|0.0003478|0.004643|0.0003236|0|0|0|1.886e-05|0|0|0.004998|AA||||||||,'.
  'T|upstream_gene_variant|MODIFIER|AP000223.42|ENSG00000260583|Transcript|ENST00000567517|antisense||||||||||'.
  'rs142513484|2407|-1||SNV|Clone_based_vega_gene||YES|||||||||||||0.0010|0.003|0.0014|0|0|0|0.004998|0|'.
  '0.0003478|0.004643|0.0003236|0|0|0|1.886e-05|0|0|0.004998|AA||||||||,T|regulatory_region_variant|MODIFIER|'.
  '||RegulatoryFeature|ENSR00001963192|TF_binding_site||||||||||rs142513484||||SNV||||||||||||||||0.0010|0.003|'.
  '0.0014|0|0|0|0.004998|0|0.0003478|0.004643|0.0003236|0|0|0|1.886e-05|0|0|0.004998|AA||||||||'.
  "\tGT\t0|0",
  'get_all_lines_by_InputBuffer - everything'
);

# custom
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 2 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

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


  $runner = get_runner({
    input_data => "21\t25585733\t.\tCATG\tTACG",
    custom => [$test_cfg->{custom_vcf}.',test,vcf,overlap'],
    output_format => 'vcf',
  });
  $of = $runner->get_OutputFactory;

  is(
    $of->get_all_lines_by_InputBuffer($runner->get_InputBuffer)->[0],
    "21\t25585733\t.\tCATG\tTACG\t.\t.\t".
    "CSQ=TACG|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1119-1122|||||||-1|||test1&del1&del2,".
    "TACG|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1030-1033|988-991|330-331|HA/RT|CATGca/CGTAca|||-1|||test1&del1&del2,".
    "TACG|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|||test1&del1&del2",
    'get_all_lines_by_InputBuffer - custom overlap'
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

# test converting to VCF from different input
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 25585733 C/- 1)]),
  dir => $test_cfg->{cache_root_dir},
  fasta => $test_cfg->{fasta},
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});

is_deeply(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585732\t21_25585733_C/-\tGC\tG\t.\t.\t".
  'CSQ=-|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,-|frameshift_variant|HIGH||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/X|Gca/ca|||-1|,-|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  "non-VCF input - deletion"
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
  "trash existing CSQ 1"
);

$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 . C T . . BAR=blah;CSQ=foo)]),
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});

is_deeply(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585733\t.\tC\tT\t.\t.\t".
  'BAR=blah;CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  "trash existing CSQ 2"
);

$of->{keep_csq} = 1;
is_deeply(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585733\t.\tC\tT\t.\t.\t".
  'BAR=blah;CSQ=foo;CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  "keep existing CSQ"
);
$of->{keep_csq} = 0;

# check we dont trash BCSQ
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 . C T . . BAR=blah;BCSQ=foo)]),
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});

is_deeply(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585733\t.\tC\tT\t.\t.\t".
  'BAR=blah;BCSQ=foo;CSQ=T|3_prime_UTR_variant|MODIFIER||ENSG00000154719|Transcript|ENST00000307301||||||1122|||||||-1|,T|missense_variant|MODERATE||ENSG00000154719|Transcript|ENST00000352957||||||1033|991|331|A/T|Gca/Aca|||-1|,T|upstream_gene_variant|MODIFIER||ENSG00000260583|Transcript|ENST00000567517||||||||||||2407|-1|',
  "dont trash BCSQ"
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


# no cons
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 25585733 C/T 1)]),
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});
$ib->buffer->[0]->{transcript_variations} = {};

is(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  "21\t25585733\t21_25585733_C/T\tC\tT\t.\t.\t.",
  "no consequences added"
);


# SV converted to VCF retains END
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(21 25585733 25585735 DEL 1)]),
  dir => $test_cfg->{cache_root_dir},
  fasta => $test_cfg->{fasta},
})->get_InputBuffer;
$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});
$ib->buffer->[0]->{transcript_variations} = {};

ok(
  $of->get_all_lines_by_InputBuffer($ib)->[0] =~ /END=25585735;CSQ/,
  "SV converted to VCF retains END"
);


## test getting stuff from input
################################

$runner = get_runner({
  input_file => $test_cfg->{test_vcf},
  dir => $test_cfg->{cache_root_dir},
  vcf => 1,
});

$of = $runner->get_OutputFactory;

is_deeply(
  [map {$of->headers->[$_]} 0..5],
  [
    '##fileformat=VCFv4.1',
    '##contig=<ID=21,assembly=GCF_000001405.26,length=46709983>',
    '##contig=<ID=22,assembly=GCF_000001405.26,length=50818468>',
    '##ALT=<ID=CNV,Description="Copy Number Polymorphism">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
  ],
  'headers - from input 1'
);


is(
  $of->headers->[-2],
  '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">',
  'headers - from input 2'
);

is(
  $of->headers->[-1],
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



## RefSeq MT transcripts have odd names - check they are not filtered out
$ib = get_runner({
  input_file => $test_cfg->create_input_file([qw(MT 12848 rs267606899 C T . . .)]),
  refseq => 1,
  use_given_ref => 1,
  offline => 1,
  dir => $test_cfg->{cache_root_dir},
})->get_InputBuffer;

$of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({config => $ib->config});


is(
  $of->get_all_lines_by_InputBuffer($ib)->[0],
  'MT	12848	rs267606899	C	T	.	.	CSQ='.
'T|downstream_gene_variant|MODIFIER||4508|Transcript|4508||||||||||||3641|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4509|Transcript|4509||||||||||||4276|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4513|Transcript|4513||||||||||||4579|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4514|Transcript|4514||||||||||||2858|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|upstream_gene_variant|MODIFIER||4519|Transcript|4519||||||||||||1899|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4537|Transcript|4537||||||||||||2444|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4538|Transcript|4538||||||||||||711|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4539|Transcript|4539||||||||||||2082|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|missense_variant|MODERATE||4540|Transcript|4540||||||512|512|171|A/V|gCa/gTa|||1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4541|Transcript|4541||||||||||||1301|-1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4556|Transcript|4556||||||||||||1826|-1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4563|Transcript|4563||||||||||||2790|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4564|Transcript|4564||||||||||||642|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4566|Transcript|4566||||||||||||4484|1||rseq_mrna_nonmatch&rseq_no_comparison,'.
'T|downstream_gene_variant|MODIFIER||4568|Transcript|4568||||||||||||512|1||rseq_mrna_nonmatch&rseq_no_comparison,T|downstream_gene_variant|MODIFIER||4571|Transcript|4571||||||||||||3108|-1||rseq_mrna_nonmatch&rseq_no_comparison,T|downstream_gene_variant|MODIFIER||4573|Transcript|4573||||||||||||2379|1||rseq_mrna_nonmatch&rseq_no_comparison,T|downstream_gene_variant|MODIFIER||4575|Transcript|4575||||||||||||583|1||rseq_mrna_nonmatch&rseq_no_comparison,T|upstream_gene_variant|MODIFIER||4576|Transcript|4576||||||||||||3040|1||rseq_mrna_nonmatch&rseq_no_comparison',
  "RefSeq MT transcripts are returned"
);


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
