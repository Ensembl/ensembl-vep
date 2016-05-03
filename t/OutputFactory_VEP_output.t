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

is_deeply(
  $of->headers(),
  [
    '## ENSEMBL VARIANT EFFECT PREDICTOR v1',
    '## Output produced at test',
    '## Using API version 1, DB version 1',
    '## Extra column keys:',
    '## IMPACT : Subjective impact classification of consequence type',
    '## DISTANCE : Shortest distance from variant to transcript',
    '## STRAND : Strand of the feature (1/-1)',
    '## FLAGS : Transcript quality flags'
  ],
  'headers'
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

my $ib = get_annotated_buffer({input_file => $test_cfg->{test_vcf}});

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
    IMPACT=MODIFIER;STRAND=-1
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
    IMPACT=MODERATE;STRAND=-1;FLAGS=cds_start_NF
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
  'CANONICAL=YES;TSL=5;APPRIS=A2;CCDS=CCDS33522.1;ENSP=ENSP00000305682;'.
  'SWISSPROT=Q9NYK5;UNIPARC=UPI00001AEAC0;EXON=11/11;'.
  'HGVSc=ENST00000307301.11:c.*18G>A;GMAF=T:0.0010;AFR_MAF=T:0.0030;'.
  'AMR_MAF=T:0.0014;EAS_MAF=T:0.0000;EUR_MAF=T:0.0000;SAS_MAF=T:0.0000;'.
  'AA_MAF=T:0.005;EA_MAF=T:0;'.
  'ExAC_MAF=T:4.119e-04;ExAC_Adj_MAF=T:0.0004133;ExAC_AFR_MAF=T:0.004681;'.
  'ExAC_AMR_MAF=T:0.000173;ExAC_EAS_MAF=T:0;ExAC_FIN_MAF=T:0;'.
  'ExAC_NFE_MAF=T:0;ExAC_OTH_MAF=T:0;ExAC_SAS_MAF=T:0',
  'get_all_lines_by_InputBuffer - everything'
);

done_testing();

sub get_annotated_buffer {
  my $tmp_cfg = shift;

  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir}.'/sereal',
    %$tmp_cfg,
  });

  $runner->init;

  my $ib = $runner->get_InputBuffer;
  $ib->next();
  $_->annotate_InputBuffer($ib) for @{$runner->get_all_AnnotationSources};
  $ib->finish_annotation();

  return $ib;
}

