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
use_ok('Bio::EnsEMBL::VEP::OutputFactory::Tab');

use_ok('Bio::EnsEMBL::VEP::Config');
use_ok('Bio::EnsEMBL::VEP::Runner');
my $cfg = Bio::EnsEMBL::VEP::Config->new();

my $of = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({config => $cfg, header_info => $test_cfg->{header_info}});

is(ref($of), 'Bio::EnsEMBL::VEP::OutputFactory::Tab', 'check class');



## METHOD TESTS
###############

is_deeply(
  $of->fields,
  [qw(
    Uploaded_variation
    Location
    Allele
    Gene
    Feature
    Feature_type
    Consequence
    cDNA_position
    CDS_position
    Protein_position
    Amino_acids
    Codons
    Existing_variation
    IMPACT
    DISTANCE
    STRAND
    FLAGS
    custom_test
  )],
  'fields'
);

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
    '## IMPACT : Subjective impact classification of consequence type',
    '## DISTANCE : Shortest distance from variant to transcript',
    '## STRAND : Strand of the feature (1/-1)',
    '## FLAGS : Transcript quality flags',
    '## custom_test : test.vcf.gz (overlap)',
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tcustom_test"
  ],
  'headers'
);

$of = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({
  config => $cfg,
  header_info => $test_cfg->{header_info}
});
$of->param('merged', 1);
$of->param('custom', 1);
is(scalar (grep {$_ eq 'SOURCE'} @{$of->fields}), 1, '--merged and --custom dont duplicate SOURCE header');
$of->param('merged', 0);
$of->param('custom', 0);

my $runner = get_annotated_buffer_runner({
  input_file => $test_cfg->{test_vcf},
  plugin => ['TestPlugin'],
  tab => 1,
  quiet => 1,
  show_ref_allele => 1 # Include reference allele in output (and header)
});
is(
  $runner->get_OutputFactory->headers->[-2].$runner->get_OutputFactory->headers->[-1],
  "## test : header".
  "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tREF_ALLELE\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\ttest",
  'headers - plugin'
);

$of = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({
  config => $cfg,
  header_info => $test_cfg->{header_info}
});
is($of->output_hash_to_line({}), '-'.("\t\-" x 17), 'output_hash_to_line - empty');

is(
  $of->output_hash_to_line({
    Uploaded_variation => 0,
  }),
  '0'.("\t\-" x 17),
  'output_hash_to_line - test 0'
);

# Include reference allele in output
my $ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  'show_ref_allele' => 1
});
$of = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({config => $ib->config});

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
    C
    MODIFIER
    -
    -1
    -
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
    C
    MODERATE
    -
    -1
    cds_start_NF
  )),
  'get_all_lines_by_InputBuffer - check last'
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
    output_format => 'tab',
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
      MODIFIER
      -
      -1
      -
      -
      test1
      BAR
    )),
    'get_all_lines_by_InputBuffer - custom'
  );
}


$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  everything => 1,
  dir => $test_cfg->{cache_root_dir},
});
$of = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({config => $ib->config});
@lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(
  $lines[0],
  "rs142513484\t21:25585733\tT\tENSG00000154719\tENST00000307301\tTranscript\t3_prime_UTR_variant\t1122\t".
  "-\t-\t-\t-\trs142513484\tMODIFIER\t-\t-1\t-\tSNV\tMRPL39\tHGNC\tHGNC:14027\tprotein_coding\tYES\t-\t-\t5\t-\t".
  "CCDS33522.1\tENSP00000305682\tQ9NYK5\t-\tUPI00001AEAC0\t-\t-\t-\t-\t11/11\t-\t-\t-\tENST00000307301.11:c.*18G>A\t".
  "-\t-\t0.0010\t0.003\t0.0014\t0\t0\t0\t0.004998\t0\t0.0003478\t0.004643\t0.0003236\t0\t0\t0\t1.886e-05\t0\t0\t".
  "0.004998\tAA\t-\t-\t-\t-\t-\t-\t-\t-\t-",
  'get_all_lines_by_InputBuffer - everything'
);


$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  everything => 1,
  dir => $test_cfg->{cache_root_dir},
  fields => 'Location,HGVSc'
});
$of = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({config => $ib->config});
@lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(
  $lines[0],
  "21:25585733\tENST00000307301.11:c.*18G>A",
  'get_all_lines_by_InputBuffer - configure fields'
);

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

