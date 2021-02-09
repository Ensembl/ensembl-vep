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

my $cfg_hash = $test_cfg->base_testing_cfg;

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::OutputFactory::VEP_output');

use_ok('Bio::EnsEMBL::VEP::Config');
use_ok('Bio::EnsEMBL::VEP::Runner');
my $cfg = Bio::EnsEMBL::VEP::Config->new();

my $of = Bio::EnsEMBL::VEP::OutputFactory::VEP_output->new({config => $cfg, header_info => $test_cfg->{header_info}});

is(ref($of), 'Bio::EnsEMBL::VEP::OutputFactory::VEP_output', 'check class');



## METHOD TESTS
###############

is_deeply($of->fields, [qw(IMPACT DISTANCE STRAND FLAGS)], 'fields');

is_deeply($of->field_order, {
  'IMPACT' => 0,
  'DISTANCE' => 1,
  'STRAND' => 2,
  'FLAGS' => 3,
}, 'field_order');
delete($of->{field_order});
delete($of->{fields});

$of->param('sift', 'b');
is_deeply($of->field_order, {
  'IMPACT' => 0,
  'DISTANCE' => 1,
  'STRAND' => 2,
  'FLAGS' => 3,
  'SIFT' => 4,
}, 'field_order - test add flag');
$of->param('sift', 0);
delete($of->{field_order});
delete($of->{fields});

$of->param('merged', 1);
$of->param('custom', 1);
is(scalar (grep {$_ eq 'SOURCE'} @{$of->fields}), 1, '--merged and --custom dont duplicate SOURCE header');
$of->param('merged', 0);
$of->param('custom', 0);
delete($of->{field_order});
delete($of->{fields});

is_deeply(
  $of->headers(),
  [
    '## ENSEMBL VARIANT EFFECT PREDICTOR v1',
    '## Output produced at test',
    '## Using API version 1, DB version 1',
    '## Column descriptions:',
    '## Uploaded_variation : Identifier of uploaded variant',
    '## Location : Location of variant in standard coordinate format (chr:start or chr:start-end)',
    '## Allele : The variant allele used to calculate the consequence',
    '## Gene : Stable ID of affected gene',
    '## Feature : Stable ID of feature',
    '## Feature_type : Type of feature - Transcript, RegulatoryFeature or MotifFeature',
    '## Consequence : Consequence type',
    '## cDNA_position : Relative position of base pair in cDNA sequence',
    '## CDS_position : Relative position of base pair in coding sequence',
    '## Protein_position : Relative position of amino acid in protein',
    '## Amino_acids : Reference and variant amino acids',
    '## Codons : Reference and variant codon sequence',
    '## Existing_variation : Identifier(s) of co-located known variants',
    '## Extra column keys:',
    '## IMPACT : Subjective impact classification of consequence type',
    '## DISTANCE : Shortest distance from variant to transcript',
    '## STRAND : Strand of the feature (1/-1)',
    '## FLAGS : Transcript quality flags',
    '## custom_test : test.vcf.gz (overlap)',
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra"
  ],
  'headers'
);

my $runner = get_annotated_buffer_runner({
  input_file => $test_cfg->{test_vcf},
  plugin => ['TestPlugin'],
  quiet => 1,
});
is(
  $runner->get_OutputFactory->headers->[-2].$runner->get_OutputFactory->headers->[-1],
  "## test : header".
  "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra",
  'headers - plugin'
);


is($of->output_hash_to_line({}), '-'.("\t\-" x 13), 'output_hash_to_line - empty');

is(
  $of->output_hash_to_line({
    Uploaded_variation => 0,
  }),
  '0'.("\t\-" x 13),
  'output_hash_to_line - test 0'
);

is(
  $of->output_hash_to_line({
    Existing_variation => 'rs123',
    Foo => 'bar',
  }),
  '-'.("\t\-" x 11)."\trs123\tFoo\=bar",
  'output_hash_to_line - test extra 1'
);

is(
  $of->output_hash_to_line({
    Existing_variation => 'rs123',
    Foo => 'bar',
    IMPACT => 'HIGH'
  }),
  '-'.("\t\-" x 11)."\trs123\tIMPACT\=HIGH;Foo\=bar",
  'output_hash_to_line - test extra 2'
);

# Include reference allele in output
my $ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  show_ref_allele => 1
});
$of = Bio::EnsEMBL::VEP::OutputFactory::VEP_output->new({config => $ib->config});

my @lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(scalar @lines, 744, 'get_all_lines_by_InputBuffer - count');

is(
  $lines[0],
  join("\t", qw(
    rs142513484
    21:25585733
    T
    ENSG00000154719
    ENST00000307301
    Transcript
    3_prime_UTR_variant
    1122
    - - - - -
    REF_ALLELE=C;IMPACT=MODIFIER;STRAND=-1
  )),
  'get_all_lines_by_InputBuffer - check first'
);

is(
  $lines[-1],
  join("\t", qw(
    rs141331202
    21:25982445
    T
    ENSG00000142192
    ENST00000448850
    Transcript
    missense_variant
    830
    832
    278
    V/I
    Gtt/Att
    -
    REF_ALLELE=C;IMPACT=MODERATE;STRAND=-1;FLAGS=cds_start_NF
  )),
  'get_all_lines_by_InputBuffer - check last'
);


$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  everything => 1,
  dir => $test_cfg->{cache_root_dir},
});
$of = Bio::EnsEMBL::VEP::OutputFactory::VEP_output->new({config => $ib->config});
@lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(
  (split("\t", $lines[0]))[-1],
  'IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=MRPL39;'.
  'SYMBOL_SOURCE=HGNC;HGNC_ID=HGNC:14027;BIOTYPE=protein_coding;'.
  'CANONICAL=YES;TSL=5;CCDS=CCDS33522.1;ENSP=ENSP00000305682;'.
  'SWISSPROT=Q9NYK5;UNIPARC=UPI00001AEAC0;EXON=11/11;'.
  'HGVSc=ENST00000307301.11:c.*18G>A;'.
  'AF=0.0010;AFR_AF=0.003;AMR_AF=0.0014;EAS_AF=0;EUR_AF=0;SAS_AF=0;AA_AF=0.004998;EA_AF=0;'.
  'gnomAD_AF=0.0003478;gnomAD_AFR_AF=0.004643;gnomAD_AMR_AF=0.0003236;gnomAD_ASJ_AF=0;'.
  'gnomAD_EAS_AF=0;gnomAD_FIN_AF=0;gnomAD_NFE_AF=1.886e-05;gnomAD_OTH_AF=0;gnomAD_SAS_AF=0;MAX_AF=0.004998;MAX_AF_POPS=AA',
  'get_all_lines_by_InputBuffer - everything'
);


# custom
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 1 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

  $runner = get_annotated_buffer_runner({
    input_file => $test_cfg->{test_vcf},
    custom => [$test_cfg->{custom_vcf}.',test,vcf,exact,,FOO'],
    output_format => 'vep',
  });
  $of = $runner->get_OutputFactory;

  @lines = @{$of->get_all_lines_by_InputBuffer($runner->get_InputBuffer)};

  is(
    $lines[0],
    join("\t", qw(
      rs142513484
      21:25585733
      T
      ENSG00000154719
      ENST00000307301
      Transcript
      3_prime_UTR_variant
      1122
      - - - - -
      IMPACT=MODIFIER;STRAND=-1;test=test1;test_FILTER=PASS;test_FOO=BAR
    )),
    'get_all_lines_by_InputBuffer - custom'
  );
}

done_testing();

sub get_annotated_buffer {
  my $tmp_cfg = shift;

  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir},
    %$tmp_cfg,
  });

  $runner->init;

  my $ib = $runner->get_InputBuffer;
  $ib->next();
  $_->annotate_InputBuffer($ib) for @{$runner->get_all_AnnotationSources};
  $ib->finish_annotation();

  return $ib;
}

sub get_annotated_buffer_runner {
  my $tmp_cfg = shift;

  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir},
    %$tmp_cfg,
  });

  $runner->init;

  my $ib = $runner->get_InputBuffer;
  $ib->next();
  $_->annotate_InputBuffer($ib) for @{$runner->get_all_AnnotationSources};
  $ib->finish_annotation();

  return $runner;
}


