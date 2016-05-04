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
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS"
  ],
  'headers'
);


is($of->output_hash_to_line({}), '-'.("\t\-" x 16), 'output_hash_to_line - empty');

is(
  $of->output_hash_to_line({
    Uploaded_variation => 0,
  }),
  '0'.("\t\-" x 16),
  'output_hash_to_line - test 0'
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
    MODERATE
    -
    -1
    cds_start_NF
  )),
  'get_all_lines_by_InputBuffer - check last'
);


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
  "-\t-\t-\t-\trs142513484\tMODIFIER\t-\t-1\t-\tSNV\tMRPL39\tHGNC\tHGNC:14027\tprotein_coding\tYES\t5\tA2\t".
  "CCDS33522.1\tENSP00000305682\tQ9NYK5\t-\tUPI00001AEAC0\t-\t-\t-\t11/11\t-\t-\tENST00000307301.11:c.*18G>A\t".
  "-\t-\tT:0.0010\tT:0.0030\tT:0.0014\tT:0.0000\tT:0.0000\tT:0.0000\tT:0.005\tT:0\tT:4.119e-04\tT:0.0004133\t".
  "T:0.004681\tT:0.000173\tT:0\tT:0\tT:0\tT:0\tT:0\t-\t-\t-\t-\t-\t-\t-\t-",
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

