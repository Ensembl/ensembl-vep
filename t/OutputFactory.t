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
use_ok('Bio::EnsEMBL::VEP::OutputFactory');


## WE NEED A RUNNER
###################

# use test
use_ok('Bio::EnsEMBL::VEP::Runner');

my $runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  input_file => $test_cfg->{test_vcf},
  check_existing => 1,
  dir => $test_cfg->{cache_root_dir},
});
ok($runner, 'new is defined');

is(ref($runner), 'Bio::EnsEMBL::VEP::Runner', 'check class');

ok($runner->init, 'init');


my $of = Bio::EnsEMBL::VEP::OutputFactory->new({config => $runner->config});
ok($of, 'new is defined');
is(ref($of), 'Bio::EnsEMBL::VEP::OutputFactory', 'check class');

my $ib = $runner->get_InputBuffer;
$ib->next();
$_->annotate_InputBuffer($ib) for @{$runner->get_all_AnnotationSources};
$ib->finish_annotation();

my $vf = $ib->buffer->[0];
my $exp = {
  'Allele' => 'T',
  'Existing_variation' => [
    'rs142513484'
  ]
};

is_deeply(
  $of->add_colocated_variant_info($vf, {Allele => 'T'}),
  $exp,
  'add_colocated_variant_info',
);

# frequency data from --check_frequency
$ib->buffer->[0]->{_freq_check_freqs} = {
  '1KG_ALL' => {
    A => 0.1
  }
};
is_deeply(
  $of->add_colocated_variant_info($ib->buffer->[0], {Allele => 'T'}),
  {
    'Allele' => 'T',
    'Existing_variation' => [
      'rs142513484',
    ],
    'FREQS' => [
      '1KG_ALL:A:0.1',
    ],
  },
  'add_colocated_variant_info - _freq_check_freqs',
);
delete($ib->buffer->[0]->{_freq_check_freqs});


# phenotype and clinsig come back without user param input
$ib = get_annotated_buffer({
  check_existing => 1,
  input_file => $test_cfg->create_input_file([qw(21 25891796 . C T . . .)])
});

is_deeply(
  $of->add_colocated_variant_info($ib->buffer->[0], {Allele => 'T'}),
  {
    'Allele' => 'T',
    'PHENO' => [
      1,
      1
    ],
    'CLIN_SIG' => [
      'uncertain_significance',
      'not_provided',
      'pathogenic'
    ],
    'Existing_variation' => [
      'rs63750066',
      'CM930033'
    ]
  },
  'add_colocated_variant_info - pheno, clin_sig, multiple',
);



## pubmed
#########

$ib = get_annotated_buffer({
  check_existing => 1,
  pubmed => 1,
  input_file => $test_cfg->create_input_file([qw(21 25272769 . C T . . .)])
});

$of->{pubmed} = 1;
is_deeply(
  $of->add_colocated_variant_info($ib->buffer->[0], {Allele => 'T'}),
  {
    'Allele' => 'T',
    'Existing_variation' => [
      'rs9977253',
    ],
    'PHENO' => [
      1
    ],
    'PUBMED' => [
      '20708005'
    ],
  },
  'add_colocated_variant_info - pubmed',
);
$of->{pubmed} = 0;


## somatic
##########

$ib = get_annotated_buffer({
  check_existing => 1,
  input_file => $test_cfg->create_input_file([qw(21 25891785 . G A . . .)])
});

is_deeply(
  $of->add_colocated_variant_info($ib->buffer->[0], {Allele => 'A'}),
  {
    'Allele' => 'A',
    'Existing_variation' => [
      'rs145564988',
      'COSM1029633',
    ],
    'SOMATIC' => [
      0, 1,
    ],
    'PHENO' => [
      0, 1,
    ],
  },
  'add_colocated_variant_info - somatic',
);



## VariationFeature_to_output_hash
##################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
});

is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'rs142513484',
    'Location' => '21:25585733'
  },
  'VariationFeature_to_output_hash'
);

$of->{variant_class} = 1;
is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'rs142513484',
    'Location' => '21:25585733',
    'VARIANT_CLASS' => 'SNV',
  },
  'VariationFeature_to_output_hash - variant_class'
);
$of->{variant_class} = 0;

$ib->buffer->[0]->{overlapping_svs} = {'sv1' => 1, 'sv2' => 2};
is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'rs142513484',
    'Location' => '21:25585733',
    'SV' => ['sv1', 'sv2'],
  },
  'VariationFeature_to_output_hash - overlapping_svs'
);


no warnings 'qw';
$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([
    ['##fileformat=VCFv4.1'],
    [qw(#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT dave barry jeff)],
    [qw(21 25607429 indtest A G . . . GT 0|1 1/1 0/0)],
  ]),
  individual => 'dave',
});


$of->{individual} = ['dave'];
is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'indtest',
    'Location' => '21:25607429',
    'IND' => 'dave',
    'ZYG' => 'HET',
  },
  'VariationFeature_to_output_hash - individual'
);
delete($of->{individual});


$ib->buffer->[0]->{_custom_annotations} = {
  custom1 => [ {name => 'foo'} ],
  custom2 => [ {name => 'bar'}, {name => 'car'} ],
  custom3 => [ {allele => 'A', name => 'shoe'} ], # this one will be ignored as it has allele key
};
is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'indtest',
    'Location' => '21:25607429',
    'IND' => 'dave',
    'ZYG' => 'HET',
    'custom1' => ['foo'],
    'custom2' => ['bar', 'car'],
  },
  'VariationFeature_to_output_hash - custom annotations'
);
delete($ib->buffer->[0]->{_custom_annotations});

$ib->buffer->[0]->{nearest} = ['foo', 'bar'];
is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'indtest',
    'Location' => '21:25607429',
    'IND' => 'dave',
    'ZYG' => 'HET',
    'NEAREST' => ['foo', 'bar'],
  },
  'VariationFeature_to_output_hash - nearest'
);
delete($ib->buffer->[0]->{nearest});

$of->{ambiguity} = 1;
is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'indtest',
    'Location' => '21:25607429',
    'AMBIGUITY' => 'R',
    'IND' => 'dave',
    'ZYG' => 'HET',
  },
  'VariationFeature_to_output_hash - ambiguity'
);
$of->{ambiguity} = 0;


## pick_worst_VariationFeatureOverlapAllele
###########################################

is($of->pick_worst_VariationFeatureOverlapAllele([]), undef, 'pick_worst_VariationFeatureOverlapAllele - empty list');

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  regulatory => 1,
});
my @vfoas =
  map {@{$_->get_all_alternate_VariationFeatureOverlapAlleles}}
  @{$ib->buffer->[0]->get_all_VariationFeatureOverlaps};

is(
  $of->pick_worst_VariationFeatureOverlapAllele(\@vfoas)->feature->stable_id,
  'ENST00000307301',
  'pick_worst_VariationFeatureOverlapAllele - default'
);

my $orig_order = $of->{pick_order};

$of->{pick_order} = ['rank'];
is(
  $of->pick_worst_VariationFeatureOverlapAllele(\@vfoas)->feature->stable_id,
  'ENST00000352957',
  'pick_worst_VariationFeatureOverlapAllele - rank'
);

$of->{pick_order} = ['appris'];
is(
  $of->pick_worst_VariationFeatureOverlapAllele(\@vfoas)->feature->stable_id,
  'ENST00000352957',
  'pick_worst_VariationFeatureOverlapAllele - appris'
);

$of->{pick_order} = ['canonical','biotype'];
is(
  $of->pick_worst_VariationFeatureOverlapAllele(\@vfoas)->feature->stable_id,
  'ENST00000307301',
  'pick_worst_VariationFeatureOverlapAllele - canonical,biotype'
);

$of->{pick_order} = ['canonical'];
is(
  $of->pick_worst_VariationFeatureOverlapAllele(\@vfoas)->feature->stable_id,
  'ENST00000307301',
  'pick_worst_VariationFeatureOverlapAllele - canonical - bail out'
);

$of->{pick_order} = $orig_order;


## pick_VariationFeatureOverlapAllele_per_gene
##############################################

is_deeply(
  [sort map {$_->feature->stable_id} @{$of->pick_VariationFeatureOverlapAllele_per_gene(\@vfoas)}],
  ['ENST00000307301', 'ENST00000567517'],
  'pick_VariationFeatureOverlapAllele_per_gene'
);

$of->{pick_order} = ['rank'];
is_deeply(
  [sort map {$_->feature->stable_id} @{$of->pick_VariationFeatureOverlapAllele_per_gene(\@vfoas)}],
  ['ENST00000352957', 'ENST00000567517'],
  'pick_VariationFeatureOverlapAllele_per_gene - change pick_order'
);
$of->{pick_order} = $orig_order;


## filter_VariationFeatureOverlapAlleles
########################################

is_deeply($of->filter_VariationFeatureOverlapAlleles([]), [], 'filter_VariationFeatureOverlapAlleles - empty arrayref');

is(scalar @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}, scalar @vfoas, 'filter_VariationFeatureOverlapAlleles - no filter');

$of->{pick} = 1;
is_deeply(
  [sort map {$_->feature->stable_id} @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}],
  ['ENST00000307301'],
  'filter_VariationFeatureOverlapAlleles - pick'
);
$of->{pick} = 0;

$of->{per_gene} = 1;
is_deeply(
  [sort map {$_->feature->stable_id} @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}],
  ['ENST00000307301', 'ENST00000567517'],
  'filter_VariationFeatureOverlapAlleles - per_gene'
);
$of->{per_gene} = 0;

$of->{flag_pick} = 1;
is(
  scalar @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)},
  scalar @vfoas,
  'filter_VariationFeatureOverlapAlleles - flag_pick count'
);
is_deeply(
  [
    sort
    map {$_->feature->stable_id}
    grep {$_->{PICK}}
    @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}
  ],
  ['ENST00000307301'],
  'filter_VariationFeatureOverlapAlleles - flag_pick check'
);
$of->{flag_pick} = 0;

# per allele tests
$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25585733 rs142513484 C G,T . . .)]),
  regulatory => 1,
});
@vfoas =
  map {@{$_->get_all_alternate_VariationFeatureOverlapAlleles}}
  @{$ib->buffer->[0]->get_all_VariationFeatureOverlaps};

$of->{pick_allele} = 1;
is_deeply(
  [
    sort
    map {$_->variation_feature_seq.':'.$_->feature->stable_id}
    @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}
  ],
  ['G:ENST00000307301', 'T:ENST00000307301'],
  'filter_VariationFeatureOverlapAlleles - pick_allele'
);
$of->{pick_allele} = 0;

$of->{flag_pick_allele} = 1;
is(
  scalar @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)},
  scalar @vfoas,
  'filter_VariationFeatureOverlapAlleles - flag_pick_allele count'
);
is_deeply(
  [
    sort
    map {$_->variation_feature_seq.':'.$_->feature->stable_id}
    grep {$_->{PICK}}
    @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}
  ],
  ['G:ENST00000307301', 'T:ENST00000307301'],
  'filter_VariationFeatureOverlapAlleles - flag_pick_allele check'
);
$of->{flag_pick_allele} = 0;

$of->{pick_allele_gene} = 1;
is_deeply(
  [
    sort
    map {$_->variation_feature_seq.':'.$_->feature->stable_id}
    @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}
  ],
  [
    'G:ENST00000307301', 'G:ENST00000567517',
    'T:ENST00000307301', 'T:ENST00000567517'
  ],
  'filter_VariationFeatureOverlapAlleles - pick_allele_gene'
);
$of->{pick_allele_gene} = 0;

$of->{flag_pick_allele_gene} = 1;
is(
  scalar @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)},
  scalar @vfoas,
  'filter_VariationFeatureOverlapAlleles - flag_pick_allele_gene count'
);
is_deeply(
  [
    sort
    map {$_->variation_feature_seq.':'.$_->feature->stable_id}
    grep {$_->{PICK}}
    @{$of->filter_VariationFeatureOverlapAlleles(\@vfoas)}
  ],
  [
    'G:ENST00000307301', 'G:ENST00000567517',
    'T:ENST00000307301', 'T:ENST00000567517'
  ],
  'filter_VariationFeatureOverlapAlleles - flag_pick_allele_gene check'
);
$of->{flag_pick_allele_gene} = 0;



## get_all_VariationFeatureOverlapAlleles
#########################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
});

is(
  scalar @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])},
  3,
  'get_all_VariationFeatureOverlapAlleles'
);

$of->{coding_only} = 1;
is(
  scalar @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])},
  1,
  'get_all_VariationFeatureOverlapAlleles - coding_only'
);
$of->{coding_only} = 0;

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25832817 . C A . . .)]),
});

is(
  scalar @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])},
  1,
  'get_all_VariationFeatureOverlapAlleles - no_intergenic off'
);

$of->{no_intergenic} = 1;
is(
  scalar @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])},
  0,
  'get_all_VariationFeatureOverlapAlleles - no_intergenic on'
);
$of->{no_intergenic} = 0;


$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  individual => 'all',
});

is(
  scalar @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])},
  0,
  'get_all_VariationFeatureOverlapAlleles - process_ref_homs off'
);

$of->{process_ref_homs} = 1;
is(
  scalar @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])},
  3,
  'get_all_VariationFeatureOverlapAlleles - process_ref_homs on'
);
$of->{process_ref_homs} = 0;



## summary_only
###############

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
});

# most_severe
$of->{most_severe} = 1;
is_deeply(
  $of->summary_only($ib->buffer->[0], {}, $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])),
  [
    {
      'Consequence' => [
        'missense_variant'
      ]
    }
  ],
  'summary_only - most_severe'
);
$of->{most_severe} = 0;

# summary
$of->{summary} = 1;
is_deeply(
  $of->summary_only($ib->buffer->[0], {}, $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])),
  [
    {
      'Consequence' => [
        'missense_variant',
        '3_prime_UTR_variant',
        'upstream_gene_variant'
      ]
    }
  ],
  'summary_only - most_severe'
);
$of->{summary} = 0;



## VariationFeatureOverlapAllele_to_output_hash
###############################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
});

my $vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];

is_deeply(
  $of->VariationFeatureOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      '3_prime_UTR_variant'
    ],
    'Allele' => 'T'
  },
  'VariationFeatureOverlapAllele_to_output_hash'
);

$of->{allele_number} = 1;
is_deeply(
  $of->VariationFeatureOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      '3_prime_UTR_variant'
    ],
    'Allele' => 'T',
    'ALLELE_NUM' => 1,
  },
  'VariationFeatureOverlapAllele_to_output_hash - allele_number'
);
$of->{allele_number} = 0;

$of->{flag_pick} = 1;
($vfoa) = grep {$_->{PICK}} @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])};
is_deeply(
  $of->VariationFeatureOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      '3_prime_UTR_variant'
    ],
    'Allele' => 'T',
    'PICK' => 1,
  },
  'VariationFeatureOverlapAllele_to_output_hash - pick'
);
$of->{flag_pick} = 0;

$of->{hgvsg} = 1;
is_deeply(
  $of->VariationFeatureOverlapAllele_to_output_hash($vfoa, {}, $ib->buffer->[0]),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      '3_prime_UTR_variant'
    ],
    'Allele' => 'T',
    'PICK' => 1,
    'HGVSg' => '21:g.25585733C>T',
  },
  'VariationFeatureOverlapAllele_to_output_hash - hgvsg'
);

# $of->chromosome_synonyms($test_cfg->{chr_synonyms});
# delete $ib->buffer->[0]->{_hgvs_genomic};

# is(
#   $of->VariationFeatureOverlapAllele_to_output_hash($vfoa, {}, $ib->buffer->[0])->{HGVSg},
#   'NC_000021.9:g.25585733C>T',
#   'VariationFeatureOverlapAllele_to_output_hash - hgvsg with synonyms'
# );

$of->{hgvsg} = 0;


$ib->buffer->[0]->{_custom_annotations} = {
  custom1 => [ {allele => 'T', name => 'foo'} ],
  custom2 => [ {name => 'bar'}, {name => 'car'} ], # this one will be ignored as it has no allele key
  custom3 => [ {allele => 'T', name => 'shoe'}, {allele => 'T', name => 'moo'} ],
};
is_deeply(
  $of->VariationFeatureOverlapAllele_to_output_hash($vfoa, {}, $ib->buffer->[0]),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      '3_prime_UTR_variant'
    ],
    'Allele' => 'T',
    'PICK' => 1,
    'custom1' => ['foo'],
    'custom3' => ['shoe', 'moo'],
  },
  'VariationFeatureOverlapAllele_to_output_hash - custom annotations'
);



## frequency tests
##################

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  check_existing => 1,
});

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];

$of->{af} = 1;

is_deeply(
  $of->add_colocated_frequency_data($vf, {Allele => 'T'}, $vf->{existing}->[0]),
  {Allele => 'T', AF => ['0.0010']},
  'add_colocated_frequency_data - af',
);

is_deeply(
  $of->add_colocated_frequency_data($vf, {Allele => 'G'}, $vf->{existing}->[0]),
  {Allele => 'G'},
  'add_colocated_frequency_data - af absent allele',
);

is_deeply(
  $of->add_colocated_frequency_data($vf, {Allele => 'C'}, $vf->{existing}->[0]),
  {Allele => 'C', AF => ['0.999']},
  'add_colocated_frequency_data - af other allele',
);

$of->{af} = 0;

# 1kg
$of->{af_1kg} = 1;

is_deeply(
  $of->add_colocated_frequency_data($vf, {Allele => 'T'}, $vf->{existing}->[0]),
  {
    Allele => 'T',
    AFR_AF => ['0.003'],
    AMR_AF => ['0.0014'],
    EAS_AF => ['0'],
    EUR_AF => ['0'],
    SAS_AF => ['0'],
  },
  'add_colocated_frequency_data - af_1kg',
);

$of->{af_1kg} = 0;


$of->{af} = 1;
$ib = get_annotated_buffer({
  check_existing => 1,
  input_file => $test_cfg->create_input_file([qw(21 25585733 25585733 G/A -)])
});
is_deeply(
  $of->add_colocated_frequency_data($ib->buffer->[0], {Allele => 'A'}, $ib->buffer->[0]->{existing}->[0]),
  {Allele => 'A', AF => ['0.0010']},
  'add_colocated_frequency_data - rev strand',
);
$of->{af} = 0;

$ib = get_annotated_buffer({
  check_existing => 1,
  input_file => $test_cfg->create_input_file([qw(21 25891796 . C T . . .)]),
  dir => $test_cfg->{exac_root_dir},
});
$of->{af_exac} = 1;

is_deeply(
  $of->add_colocated_frequency_data($ib->buffer->[0], {Allele => 'T'}, $ib->buffer->[0]->{existing}->[0]),
  {
    Allele => 'T',
    'ExAC_OTH_AF' => [
      '0.001101'
    ],
    'ExAC_Adj_AF' => [
      '5.768e-05'
    ],
    'ExAC_AFR_AF' => [
      '0'
    ],
    'ExAC_AMR_AF' => [
      '0.0003457'
    ],
    'ExAC_NFE_AF' => [
      '2.998e-05'
    ],
    'ExAC_SAS_AF' => [
      '0'
    ],
    'ExAC_FIN_AF' => [
      '0'
    ],
    'ExAC_EAS_AF' => [
      '0'
    ],
    'ExAC_AF' => [
      '5.765e-05'
    ]
  },
  'add_colocated_frequency_data - af_exac',
);
$of->{af_exac} = 0;

$ib = get_annotated_buffer({
  check_existing => 1,
  input_file => $test_cfg->create_input_file([qw(21 25891796 . C T . . .)]),
});
$of->{af_gnomad} = 1;

is_deeply(
  $of->add_colocated_frequency_data($ib->buffer->[0], {Allele => 'T'}, $ib->buffer->[0]->{existing}->[0]),
  {
    'Allele' => 'T',
    'gnomAD_NFE_AF' => [
      '2.687e-05'
    ],
    'gnomAD_AF' => [
      '9.75e-05'
    ],
    'gnomAD_OTH_AF' => [
      '0.0001823'
    ],
    'gnomAD_AFR_AF' => [
      '0'
    ],
    'gnomAD_EAS_AF' => [
      '0'
    ],
    'gnomAD_AMR_AF' => [
      '0.0005957'
    ],
    'gnomAD_SAS_AF' => [
      '0'
    ],
    'gnomAD_FIN_AF' => [
      '0'
    ],
    'gnomAD_ASJ_AF' => [
      '0'
    ]
  },
  'add_colocated_frequency_data - af_gnomad',
);
$of->{af_gnomad} = 0;

$ib = get_annotated_buffer({
  check_existing => 1,
  input_file => $test_cfg->create_input_file([qw(21 25975223 . G A . . .)])
});

$of->{af_esp} = 1;
is_deeply(
  $of->add_colocated_frequency_data($ib->buffer->[0], {Allele => 'A'}, $ib->buffer->[0]->{existing}->[0]),
  {
    'Allele' => 'A',
    'AA_AF' => [
      '0',
    ],
    'EA_AF' => [
      '0.000814',
    ],
  },
  'add_colocated_frequency_data - af_esp',
);
$of->{af_esp} = 0;


# max_af
$ib = get_annotated_buffer({
  check_existing => 1,
  input_file => $test_cfg->create_input_file([qw(21 25594005 . A G . . .)])
});

$of->{max_af} = 1;
is_deeply(
  $of->add_colocated_frequency_data($ib->buffer->[0], {Allele => 'G'}, $ib->buffer->[0]->{existing}->[0]),
  {
    'Allele' => 'G',
    'MAX_AF' => '0.8223',
    'MAX_AF_POPS' => [
      'gnomAD_ASJ',
    ],
  },
  'add_colocated_frequency_data - max_af',
);
$of->{max_af} = 0;



## BaseTranscriptVariationAllele_to_output_hash
###############################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
});

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[2];

is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa, {Consequence => ['upstream_gene_variant']}),
  {
    'STRAND' => -1,
    'Consequence' => [
      'upstream_gene_variant'
    ],
    'Feature_type' => 'Transcript',
    'Feature' => 'ENST00000567517',
    'Gene' => 'ENSG00000260583',
    'DISTANCE' => 2407,
  },
  'BaseTranscriptVariationAllele_to_output_hash'
);

$of->{transcript_version} = 1;
is(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa, {Consequence => ['upstream_gene_variant']})->{Feature},
  'ENST00000567517.1',
  'BaseTranscriptVariationAllele_to_output_hash'
);
$of->{transcript_version} = 0;

($vf) = grep {$_->{variation_name} eq 'rs199510789'} @{$ib->buffer};
($vfoa) = grep {$_->feature->stable_id eq 'ENST00000419219'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};

is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa),
  {
    'STRAND' => -1,
    'Feature_type' => 'Transcript',
    'Feature' => 'ENST00000419219',
    'Gene' => 'ENSG00000154719',
    'FLAGS' => ['cds_end_NF'],
  },
  'BaseTranscriptVariationAllele_to_output_hash - check transcript FLAGS'
);

$of->{numbers} = 1;
$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[1];
is(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{EXON},
  '10/10',
  'BaseTranscriptVariationAllele_to_output_hash - exon numbers'
);

($vf) = grep {$_->{variation_name} eq 'rs187353664'} @{$ib->buffer};
($vfoa) = grep {$_->feature->stable_id eq 'ENST00000352957'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};

is(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{INTRON},
  '9/9',
  'BaseTranscriptVariationAllele_to_output_hash - intron numbers'
);
$of->{numbers} = 0;

$of->{domains} = 1;
($vf) = grep {$_->{variation_name} eq 'rs116645811'} @{$ib->buffer};
($vfoa) = grep {$_->feature->stable_id eq 'ENST00000307301'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};
is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{DOMAINS},
  ['Low_complexity_(Seg):seg'],
  'BaseTranscriptVariationAllele_to_output_hash - domains'
);
$of->{domains} = 0;

$of->{symbol} = 1;
$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];
is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa),
  {
    'STRAND' => -1,
    'HGNC_ID' => 'HGNC:14027',
    'SYMBOL' => 'MRPL39',
    'Feature_type' => 'Transcript',
    'SYMBOL_SOURCE' => 'HGNC',
    'Gene' => 'ENSG00000154719',
    'Feature' => 'ENST00000307301'
  },
  'BaseTranscriptVariationAllele_to_output_hash - symbol'
);
$of->{symbol} = 0;

$of->{gene_phenotype} = 1;
($vf) = grep {$_->{variation_name} eq 'rs145277462'} @{$ib->buffer};

($vfoa) = grep {$_->feature->stable_id eq 'ENST00000346798'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};
is(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{GENE_PHENO},
  1,
  'BaseTranscriptVariationAllele_to_output_hash - gene_phenotype'
);
$of->{gene_phenotype} = 0;

# we can test these ones en-masse
my @flags = (
  [qw(ccds        CCDS      CCDS33522.1)],
  [qw(protein     ENSP      ENSP00000305682)],
  [qw(canonical   CANONICAL YES)],
  [qw(biotype     BIOTYPE   protein_coding)],
  [qw(tsl         TSL       5)],
);
my $method = 'BaseTranscriptVariationAllele_to_output_hash';
$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];
foreach my $flag(@flags) {
  $of->{$flag->[0]} = 1;
  is($of->$method($vfoa)->{$flag->[1]}, $flag->[2], $method.' - '.$flag->[0]);
  $of->{$flag->[0]} = 0;
}

# and these, but we expect arrayrefs
@flags = (
  [qw(xref_refseq RefSeq    NM_080794.3)],
  [qw(uniprot     SWISSPROT Q9NYK5)],
  [qw(uniprot     UNIPARC   UPI00001AEAC0)],
);
$method = 'BaseTranscriptVariationAllele_to_output_hash';
$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];

foreach my $flag(@flags) {
  $of->{$flag->[0]} = 1;
  is_deeply($of->$method($vfoa)->{$flag->[1]}, [$flag->[2]], $method.' - '.$flag->[0]);
  $of->{$flag->[0]} = 0;
}

# test miRNA
$of->{mirna} = 1;

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([
    [qw(21 25573985 25573985 C T . . .)],
    [qw(21 25573987 25573987 C CG . . .)],
    [qw(21 25573994 25573994 C T . . .)]
  ])
});

$vf = $ib->buffer->[0];
($vfoa) = grep {$_->feature->stable_id eq 'ENST00000385060'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};
is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{miRNA},
  ['miRNA_stem'],
  'BaseTranscriptVariationAllele_to_output_hash - miRNA stem'
);

$vf = $ib->buffer->[1];
($vfoa) = grep {$_->feature->stable_id eq 'ENST00000385060'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};
is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{miRNA},
  ['miRNA_stem'],
  'BaseTranscriptVariationAllele_to_output_hash - miRNA stem insertion'
);

$vf = $ib->buffer->[2];
($vfoa) = grep {$_->feature->stable_id eq 'ENST00000385060'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};
is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{miRNA},
  ['miRNA_loop'],
  'BaseTranscriptVariationAllele_to_output_hash - miRNA loop'
);

$of->{mirna} = 0;



# REFSEQ_MATCH - use refseq cache
$of->{refseq} = 1;

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  species => 'homo_sapiens',
  dir => $test_cfg->{cache_root_dir},
  refseq => 1,
  fasta => $test_cfg->{fasta},
});

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];

is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{REFSEQ_MATCH},
  [
    'rseq_mrna_nonmatch',
    'rseq_cds_mismatch'
  ],
  'BaseTranscriptVariationAllele_to_output_hash - REFSEQ_MATCH'
);
$of->{refseq} = 0;

# REFSEQ_MATCH - use merged cache
$of->{merged} = 1;

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  species => 'homo_sapiens',
  dir => $test_cfg->{cache_root_dir},
  merged => 1,
});

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[3];

is_deeply(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{REFSEQ_MATCH},
  [
    'rseq_ens_match_cds',
    'rseq_mrna_nonmatch',
    'rseq_cds_mismatch'
  ],
  'BaseTranscriptVariationAllele_to_output_hash - REFSEQ_MATCH merged'
);

is(
  $of->BaseTranscriptVariationAllele_to_output_hash($vfoa)->{SOURCE},
  'RefSeq',
  'BaseTranscriptVariationAllele_to_output_hash - merged source 1'
);

is(
  $of->BaseTranscriptVariationAllele_to_output_hash(
    $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0]
  )->{SOURCE},
  'Ensembl',
  'BaseTranscriptVariationAllele_to_output_hash - merged source 2'
);
$of->{merged} = 0;



## TranscriptVariationAllele_to_output_hash
###########################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  dir => $test_cfg->{cache_root_dir},
});

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[1];

is_deeply(
  $of->TranscriptVariationAllele_to_output_hash($vfoa, {}),
  {
    'STRAND' => -1,
    'IMPACT' => 'MODERATE',
    'Consequence' => [
      'missense_variant'
    ],
    'Feature_type' => 'Transcript',
    'Allele' => 'T',
    'CDS_position' => 991,
    'Gene' => 'ENSG00000154719',
    'cDNA_position' => 1033,
    'Protein_position' => 331,
    'Amino_acids' => 'A/T',
    'Feature' => 'ENST00000352957',
    'Codons' => 'Gca/Aca'
  },
  'TranscriptVariationAllele_to_output_hash'
);

$of->{total_length} = 1;
is_deeply(
  $of->TranscriptVariationAllele_to_output_hash($vfoa, {}),
  {
    'STRAND' => -1,
    'IMPACT' => 'MODERATE',
    'Consequence' => [
      'missense_variant'
    ],
    'Feature_type' => 'Transcript',
    'Allele' => 'T',
    'CDS_position' => '991/1017',
    'Gene' => 'ENSG00000154719',
    'cDNA_position' => '1033/1110',
    'Protein_position' => '331/338',
    'Amino_acids' => 'A/T',
    'Feature' => 'ENST00000352957',
    'Codons' => 'Gca/Aca'
  },
  'TranscriptVariationAllele_to_output_hash - total_length'
);
$of->{total_length} = 0;

$of->{hgvsc} = 1;
is(
  $of->TranscriptVariationAllele_to_output_hash($vfoa)->{HGVSc},
  'ENST00000352957.8:c.991G>A',
  'TranscriptVariationAllele_to_output_hash - HGVSc'
);
$of->{hgvsc} = 0;

$of->{hgvsp} = 1;
is(
  $of->TranscriptVariationAllele_to_output_hash($vfoa)->{HGVSp},
  'ENSP00000284967.6:p.Ala331Thr',
  'TranscriptVariationAllele_to_output_hash - HGVSp'
);
$of->{hgvsp} = 0;

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];
is_deeply(
  $of->TranscriptVariationAllele_to_output_hash($vfoa, {}),
  {
    'STRAND' => -1,
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      '3_prime_UTR_variant'
    ],
    'Feature_type' => 'Transcript',
    'Allele' => 'T',
    'Gene' => 'ENSG00000154719',
    'cDNA_position' => 1122,
    'Feature' => 'ENST00000307301'
  },
  'TranscriptVariationAllele_to_output_hash - non-coding'
);

$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  dir => $test_cfg->{cache_root_dir},
});

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[1];
$of->param('shift_length', 1);
$of->param('shift_3prime', 1);
is_deeply(
  $of->TranscriptVariationAllele_to_output_hash($vfoa, {}),
  {
    'STRAND' => -1,
    'IMPACT' => 'MODERATE',
    'Consequence' => [
      'missense_variant'
    ],
    'Feature_type' => 'Transcript',
    'Allele' => 'T',
    'CDS_position' => 991,
    'Gene' => 'ENSG00000154719',
    'cDNA_position' => 1033,
    'Protein_position' => 331,
    'Amino_acids' => 'A/T',
    'Feature' => 'ENST00000352957',
    'Codons' => 'Gca/Aca',
    'SHIFT_LENGTH' => 0,
  },
  'TranscriptVariationAllele_to_output_hash -> shift_length'
);



@flags = (
  [qw(sift     p SIFT     tolerated_low_confidence)],
  [qw(sift     s SIFT     0.17)],
  [qw(sift     b SIFT     tolerated_low_confidence\(0.17\))],
  [qw(polyphen p PolyPhen benign)],
  [qw(polyphen s PolyPhen 0.001)],
  [qw(polyphen b PolyPhen benign\(0.001\))],
);
$method = 'TranscriptVariationAllele_to_output_hash';
($vf) = grep {$_->{variation_name} eq 'rs142513484'} @{$ib->buffer};
($vfoa) = grep {$_->feature->stable_id eq 'ENST00000352957'} @{$of->get_all_VariationFeatureOverlapAlleles($vf)};

foreach my $flag(@flags) {
  $of->{$flag->[0]} = $flag->[1];
  is($of->$method($vfoa)->{$flag->[2]}, $flag->[3], $method.' - '.$flag->[0].' '.$flag->[1]);
  $of->{$flag->[0]} = 0;
}



## RegulatoryFeatureVariationAllele_to_output_hash
##################################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25734924 . C T . . .)]),
  regulatory => 1,
});

($vfoa) = grep {ref($_) =~ /Regulatory/} @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])};

is_deeply(
  $of->RegulatoryFeatureVariationAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'regulatory_region_variant'
    ],
    'Feature_type' => 'RegulatoryFeature',
    'BIOTYPE' => 'promoter',
    'Feature' => 'ENSR00000140763',
    'Allele' => 'T'
  },
  'RegulatoryFeatureVariationAllele_to_output_hash'
);

$of->{cell_type} = ['HUVEC'];
is_deeply(
  $of->RegulatoryFeatureVariationAllele_to_output_hash($vfoa)->{CELL_TYPE},
  ['HUVEC:INACTIVE'],
  'RegulatoryFeatureVariationAllele_to_output_hash - cell_type'
);
$of->{cell_type} = undef;




## MotifFeatureVariationAllele_to_output_hash
#############################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25735354 . C T . . .)]),
  regulatory => 1,
});

($vfoa) = grep {ref($_) =~ /Motif/} @{$of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])};

is_deeply(
  $of->MotifFeatureVariationAllele_to_output_hash($vfoa),
  {
    'STRAND' => 1,
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'TF_binding_site_variant'
    ],
    'MOTIF_POS' => 7,
    'Feature_type' => 'MotifFeature',
    'MOTIF_NAME' => 'ENSPFM0042',
    'Allele' => 'T',
    'TRANSCRIPTION_FACTORS' => ['CTCF'],
    'Feature' => 'ENSM00191005622',
    'HIGH_INF_POS' => 'N',
    'MOTIF_SCORE_CHANGE' => '-0.026'
  },
  'MotifFeatureVariationAllele_to_output_hash'
);

# link to cell types for new motif_feature schema has not been added yet
#$of->{cell_type} = ['MultiCell'];
#is_deeply(
#  $of->MotifFeatureVariationAllele_to_output_hash($vfoa)->{CELL_TYPE},
#  ['MultiCell:Promoter'],
#  'MotifFeatureVariationAllele_to_output_hash - cell_type'
#);
$of->{cell_type} = undef;




## IntergenicVariationAllele_to_output_hash
###########################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25832817 . C A . . .)]),
});

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];

is_deeply(
  $of->IntergenicVariationAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'intergenic_variant'
    ],
    'Allele' => 'A'
  },
  'IntergenicVariationAllele_to_output_hash'
);




##############################
##############################
#### STRUCTURAL VARIATION ####
##############################
##############################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25606614 sv_dup T . . . SVTYPE=DUP;END=25606616)]),
  regulatory => 1,
});

is_deeply(
  $of->VariationFeature_to_output_hash($ib->buffer->[0]),
  {
    'Uploaded_variation' => 'sv_dup',
    'Location' => '21:25606615-25606616'
  },
  'SV - VariationFeature_to_output_hash'
);

@vfoas =
  map {@{$_->get_all_alternate_StructuralVariationOverlapAlleles}}
  @{$ib->buffer->[0]->get_all_StructuralVariationOverlaps};

is(
  $of->pick_worst_VariationFeatureOverlapAllele(\@vfoas)->feature->stable_id,
  'ENST00000307301',
  'SV - pick_worst_VariationFeatureOverlapAllele'
);

is_deeply(
  [sort map {$_->feature->stable_id} @{$of->pick_VariationFeatureOverlapAllele_per_gene(\@vfoas)}],
  ['ENSR00000140751', 'ENST00000307301'],
  'SV - pick_VariationFeatureOverlapAllele_per_gene'
);

is(scalar @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)}, scalar @vfoas, 'SV - filter_StructuralVariationOverlapAlleles - no filter');


## filter_StructuralVariationOverlapAlleles
###########################################

# Pick
$of->{pick} = 1;
is_deeply(
  [sort map {$_->feature->stable_id} @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)}],
  ['ENST00000307301'],
  'SV - filter_StructuralVariationOverlapAlleles - pick'
);
$of->{pick} = 0;

# Pick allele
$of->{pick_allele} = 1;
is_deeply(
  [sort map {$_->feature->stable_id} @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)}],
  ['ENST00000307301'],
  'SV - filter_StructuralVariationOverlapAlleles - pic_allele'
);
$of->{pick_allele} = 0;

# Flag pick
$of->{flag_pick} = 1;
is(
  scalar @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)},
  scalar @vfoas,
  'SV - filter_StructuralVariationOverlapAlleles - flag_pick count'
);
is_deeply(
  [
    sort
    map {$_->feature->stable_id}
    grep {$_->{PICK}}
    @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)}
  ],
  ['ENST00000307301'],
  'SV - filter_StructuralVariationOverlapAlleles - flag_pick check'
);
$of->{flag_pick} = 0;

# Flag pick allele
$of->{flag_pick_allele} = 1;
is(
  scalar @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)},
  scalar @vfoas,
  'SV - filter_StructuralVariationOverlapAlleles - flag_pick_allele count'
);
is_deeply(
  [
    sort
    map {$_->feature->stable_id}
    grep {$_->{PICK}}
    @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)}
  ],
  ['ENST00000307301'],
  'SV - filter_StructuralVariationOverlapAlleles - flag_pick_allele check'
);
$of->{flag_pick_allele} = 0;

# Flag pick allele gene
$of->{flag_pick_allele_gene} = 1;
is(
  scalar @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)},
  scalar @vfoas,
  'SV - filter_StructuralVariationOverlapAlleles - flag_pick_allele_gene count'
);
is_deeply(
  [
    sort
    map {$_->feature->stable_id}
    grep {$_->{PICK}}
    @{$of->filter_StructuralVariationOverlapAlleles(\@vfoas)}
  ],
  [
    'ENSR00000140751', 'ENST00000307301',
  ],
  'SV - filter_StructuralVariationOverlapAlleles - flag_pick_allele_gene check'
);
$of->{flag_pick_allele_gene} = 0;


## get_all_StructuralVariationOverlapAlleles
############################################

is(
  scalar @{$of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])},
  4,
  'SV - get_all_StructuralVariationOverlapAlleles'
);

$of->{coding_only} = 1;
is(
  scalar @{$of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])},
  3,
  'SV - get_all_StructuralVariationOverlapAlleles - coding_only'
);
$of->{coding_only} = 0;

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25832817 . C . . . SVTYPE=DUP;END=25832818)]),
});

is(
  scalar @{$of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])},
  1,
  'SV - get_all_StructuralVariationOverlapAlleles - no_intergenic off'
);

$of->{no_intergenic} = 1;
is(
  scalar @{$of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])},
  0,
  'SV - get_all_StructuralVariationOverlapAlleles - no_intergenic on'
);
$of->{no_intergenic} = 0;


## BaseStructuralVariationOverlapAllele_to_output_hash
######################################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25585733 sv_dup T . . . SVTYPE=DUP;END=25585735)]),
});

$vfoa = $of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])->[1];

is_deeply(
  $of->BaseStructuralVariationOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'coding_sequence_variant',
      'feature_elongation'
    ],
    'Allele' => 'duplication'
  },
  'SV - BaseStructuralVariationOverlapAllele_to_output_hash'
);

$of->{allele_number} = 1;
is(
  $of->BaseStructuralVariationOverlapAllele_to_output_hash($vfoa)->{ALLELE_NUM},
  1,
  'SV - BaseStructuralVariationOverlapAllele_to_output_hash - allele_number'
);
$of->{allele_number} = 0;

## StructuralVariationOverlapAllele_to_output_hash
##################################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25585733 sv_dup T . . . SVTYPE=DUP;END=25585735)]),
});

$vfoa = $of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])->[1];

is_deeply(
  $of->StructuralVariationOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'coding_sequence_variant',
      'feature_elongation'
    ],
    'OverlapPC' => '0.01',
    'Feature_type' => 'Transcript',
    'OverlapBP' => 2,
    'Feature' => 'ENST00000352957',
    'Allele' => 'duplication'
  },
  'SV - StructuralVariationOverlapAllele_to_output_hash'
);

$of->{allele_number} = 1;
is_deeply(
  $of->StructuralVariationOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'coding_sequence_variant',
      'feature_elongation'
    ],
    'OverlapPC' => '0.01',
    'Feature_type' => 'Transcript',
    'OverlapBP' => 2,
    'Feature' => 'ENST00000352957',
    'Allele' => 'duplication',
    'ALLELE_NUM' => 1,
  },
  'SV - StructuralVariationOverlapAllele_to_output_hash'
);
$of->{allele_number} = 0;


$of->{flag_pick} = 1;
($vfoa) = grep {$_->{PICK}} @{$of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])};
is_deeply(
  $of->StructuralVariationOverlapAllele_to_output_hash($vfoa),  
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      '3_prime_UTR_variant',
      'feature_elongation'
    ],
    'OverlapPC' => '0.01',
    'Feature_type' => 'Transcript',
    'OverlapBP' => 2,
    'Feature' => 'ENST00000307301',
    'Allele' => 'duplication',
    'PICK' => 1,
  },
  'SV - StructuralVariationOverlapAllele_to_output_hash - flag_pick'
);
$of->{flag_pick} = 0;



# regulatory
$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25735354 . C . . . SVTYPE=DUP;END=25735356)]),
  regulatory => 1,
});

($vfoa) = grep {ref($_->feature) =~ /Motif/} @{$of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])};

is_deeply(
  $of->StructuralVariationOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'TF_binding_site_variant'
    ],
    'OverlapPC' => '11.76',
    'Feature_type' => 'MotifFeature',
    'OverlapBP' => 2,
    'Feature' => 'ENSM00191005622',
    'Allele' => 'duplication'
  },
  'SV - StructuralVariationOverlapAllele_to_output_hash - MotifFeature'
);

($vfoa) = grep {ref($_->feature) =~ /Regulatory/} @{$of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])};

is_deeply(
  $of->StructuralVariationOverlapAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'regulatory_region_variant'
    ],
    'OverlapPC' => '0.04',
    'Feature_type' => 'RegulatoryFeature',
    'OverlapBP' => 2,
    'Feature' => 'ENSR00000140763',
    'Allele' => 'duplication'
  },
  'SV - StructuralVariationOverlapAllele_to_output_hash - RegulatoryFeature'
);

$of->{cell_type} = ['HUVEC'];
is_deeply(
  $of->StructuralVariationOverlapAllele_to_output_hash($vfoa)->{CELL_TYPE},
  ['HUVEC:INACTIVE'],
  'SV - StructuralVariationOverlapAllele_to_output_hash - RegulatoryFeature cell_type'
);
$of->{cell_type} = undef;

# intergenic
$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25832817 . C . . . SVTYPE=DUP;END=25832818)]),
});

$vfoa = $of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])->[0];

is_deeply(
  $of->IntergenicStructuralVariationAllele_to_output_hash($vfoa),
  {
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'intergenic_variant'
    ],
    'Allele' => 'duplication'
  },
  'SV - IntergenicStructuralVariationAllele_to_output_hash'
);





## TranscriptStructuralVariationAllele_to_output_hash
#####################################################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25585733 sv_dup T . . . SVTYPE=DUP;END=25585735)]),
});

$vfoa = $of->get_all_StructuralVariationOverlapAlleles($ib->buffer->[0])->[1];

is_deeply(
  $of->TranscriptStructuralVariationAllele_to_output_hash($vfoa, {}),
  {
    'STRAND' => -1,
    'IMPACT' => 'MODIFIER',
    'Consequence' => [
      'coding_sequence_variant',
      'feature_elongation'
    ],
    'OverlapPC' => '0.01',
    'Feature_type' => 'Transcript',
    'Allele' => 'duplication',
    'CDS_position' => '989-990',
    'Gene' => 'ENSG00000154719',
    'cDNA_position' => '1031-1032',
    'Protein_position' => 330,
    'Feature' => 'ENST00000352957',
    'OverlapBP' => 2
  },
  'SV - TranscriptStructuralVariationAllele_to_output_hash'
);



## minimal restoration etc
##########################

$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25741665 . CAGAAGAAAG TAGAAGAAAG,C . . .)]),
  minimal => 1,
});

is(scalar @{$ib->buffer}, 2, 'minimal - expanded count');
is($ib->buffer->[0]->allele_string, 'C/T', 'minimal - expanded first allele string');

$of->rejoin_variants_in_InputBuffer($ib);

is(scalar @{$ib->buffer}, 1, 'minimal - rejoined count');
is($ib->buffer->[0]->allele_string, 'CAGAAGAAAG/TAGAAGAAAG/C', 'minimal - rejoined allele string');

my %by_allele =
  map {$_->{Allele} => $_}
  grep {$_->{Feature} eq 'ENST00000400075'}
  @{$of->get_all_output_hashes_by_VariationFeature($ib->buffer->[0])};

is_deeply(
  \%by_allele,
  {
    '-' => {
      'STRAND' => 1,
      'IMPACT' => 'MODERATE',
      'Consequence' => [
        'inframe_deletion',
        'splice_region_variant'
      ],
      'MINIMISED' => 1,
      'Feature_type' => 'Transcript',
      'Uploaded_variation' => '21_25741665_CAGAAGAAAG/TAGAAGAAAG/C',
      'Allele' => '-',
      'CDS_position' => '68-76',
      'Gene' => 'ENSG00000154727',
      'cDNA_position' => '289-297',
      'Protein_position' => '23-26',
      'Amino_acids' => 'KKKG/S',
      'Feature' => 'ENST00000400075',
      'Codons' => 'aAGAAGAAAGgc/agc',
      'Location' => '21:25741665-25741674',
      'SHIFT_LENGTH' => 0,
    },
    'T' => {
      'STRAND' => 1,
      'IMPACT' => 'MODERATE',
      'Consequence' => [
        'missense_variant'
      ],
      'MINIMISED' => 1,
      'Feature_type' => 'Transcript',
      'Uploaded_variation' => '21_25741665_CAGAAGAAAG/TAGAAGAAAG/C',
      'Allele' => 'T',
      'CDS_position' => 67,
      'Gene' => 'ENSG00000154727',
      'cDNA_position' => 288,
      'Protein_position' => 23,
      'Amino_acids' => 'P/S',
      'Feature' => 'ENST00000400075',
      'Codons' => 'Cca/Tca',
      'Location' => '21:25741665-25741674',
      'SHIFT_LENGTH' => 0,
    }
  },
  'minimal - output hashes'
);


# test intergenic
$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([
    [qw(21 25832817 . G GC,C . . .)],
    [qw(21 25832815 . G C . . .)]
  ]),
  minimal => 1
});

is(scalar @{$ib->buffer}, 3, 'minimal - intergenic - expanded count');

$of->rejoin_variants_in_InputBuffer($ib);

is(scalar @{$ib->buffer}, 2, 'minimal - intergenic - rejoined count');

foreach my $vf (@{$ib->buffer}) {
  my $vfoas = $of->get_all_VariationFeatureOverlapAlleles($vf);
  foreach my $vfoa (@{$vfoas}) {
    is(ref($vfoa->base_variation_feature_overlap), 'Bio::EnsEMBL::Variation::IntergenicVariation', 'base_variation_feature_overlap is set after rejoin');
  }
}

is_deeply(
  [map {$_->display_consequence} @{$ib->buffer}],
  ['intergenic_variant', 'intergenic_variant'],
  'minimal - intergenic - display_consequence check'
);


# test case where minimal resolves two ALTs to the same thing
# allele_num should keep track of them
$ib = get_annotated_buffer({
  input_file => $test_cfg->create_input_file([qw(21 25606454 test G GC,C . . .)]),
  minimal => 1
});

$of->rejoin_variants_in_InputBuffer($ib);

$of->{allele_number} = 1;

%by_allele =
  map {$_->{Consequence}->[0] => $_}
  grep {$_->{Feature} eq 'ENST00000419219'}
  @{$of->get_all_output_hashes_by_VariationFeature($ib->buffer->[0])};

is(scalar keys %{{map {$_->{Allele} => 1} values %by_allele}}, 1, "minimal - allele num where two alts resolve to same allele (check allele)");
is($by_allele{frameshift_variant}->{ALLELE_NUM}, 1, "minimal - allele num where two alts resolve to same allele (frameshift)");
is($by_allele{missense_variant}->{ALLELE_NUM}, 2, "minimal - allele num where two alts resolve to same allele (missense)");

$of->{allele_number} = 0;


## custom headers
#################

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 1 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

  $runner = get_annotated_buffer_runner({
    input_file => $test_cfg->create_input_file([qw(21 25606454 test G C . . .)]),
    custom => [$test_cfg->{custom_vcf}.',test,vcf'],
    quiet => 1,
    warning_file => 'STDERR',
  });

  $of = $runner->get_OutputFactory();
  $ib = $runner->get_InputBuffer();

  is_deeply(
    $of->get_custom_headers,
    [
      [
        'test',
        $test_cfg->{custom_vcf}.' (overlap)'
      ]
    ],
    'get_custom_headers'
  );
}


## plugins
##########

$runner = get_annotated_buffer_runner({
  input_file => $test_cfg->create_input_file([qw(21 25606454 test G C . . .)]),
  plugin => ['TestPlugin'],
  quiet => 1,
  warning_file => 'STDERR',
});

$of = $runner->get_OutputFactory();
$ib = $runner->get_InputBuffer();

is_deeply(
  $of->get_plugin_headers,
  [['test', 'header']],
  'get_plugin_headers'
);

is_deeply(
  $of->run_plugins($of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0], {}),
  {test => 'Hello'},
  'run_plugins'
);

is_deeply(
  $of->run_plugins($of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0], 'skip'),
  undef,
  'run_plugins - skip'
);


# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n";

my $tmp;
close STDERR;
open STDERR, '>', \$tmp;

is_deeply(
  $of->run_plugins($of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0], 'not_hash'),
  'not_hash',
  'run_plugins - return is not hash'
);
ok($tmp =~ /did not return a hashref/, 'run_plugins - return is not hash message');

$runner = get_annotated_buffer_runner({
  input_file => $test_cfg->create_input_file([qw(21 25606454 test G C . . .)]),
  plugin => ['TestPluginRunFails'],
  quiet => 1,
  warning_file => 'STDERR',
});

$of = $runner->get_OutputFactory();
$ib = $runner->get_InputBuffer();

is_deeply(
  $of->run_plugins($of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0], {}),
  {},
  'run_plugins - new fails'
);

ok($tmp =~ /went wrong/, 'run_plugins - new fails message');

# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";


## HGVS shifting
################

$of->{hgvs}  = 1;
$of->{hgvsg} = 1;
my $input_file_example = $test_cfg->create_input_file([qw(21 25592985 hgvsins A ATAAA . . .)]);

# Shifting ON
$ib = get_annotated_buffer({
  input_file => $input_file_example,
  shift_hgvs => 1,
  shift_length => 1, 
},1);

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];
is_deeply(
  $of->VariationFeatureOverlapAllele_to_output_hash($vfoa, {}, $ib->buffer->[0]),
  {
    "IMPACT" => "MODIFIER",
    "Consequence" => [
      "intron_variant"
    ],
    "HGVSg" => "21:g.25592986_25592989dup",
    "Allele" => "TAAA"
  },
  'HGVS 3prime shifting - ON'
);

# Shifting OFF
$ib = get_annotated_buffer({
  input_file => $input_file_example,
  shift_hgvs => 0
},1);

$vfoa = $of->get_all_VariationFeatureOverlapAlleles($ib->buffer->[0])->[0];
is_deeply(
  $of->VariationFeatureOverlapAllele_to_output_hash($vfoa, {}, $ib->buffer->[0]),
  {
    "IMPACT" => "MODIFIER",
    "Consequence" => [
      "intron_variant"
    ],
    "HGVSg" => "21:g.25592982_25592985dup",
    "Allele" => "TAAA"
  },
  'HGVS 3prime shifting - OFF'
);

$of->{hgvs}  = 0;
$of->{hgvsg} = 0;

## Check reset_shifted_positions

$vfoa->_return_3prime;
$vfoa->transcript_variation->cds_start(20);
my $new_cds_start = $vfoa->transcript_variation->{cds_start};
$of->reset_shifted_positions($vfoa->variation_feature);
ok(defined($new_cds_start) && !defined($vfoa->transcript_variation->{cds_start}), 'reset_shifted_positions');



# done
done_testing();

sub get_annotated_buffer {
  my $tmp_cfg = shift;
  my $set_package_variables = shift;

  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir},
    %$tmp_cfg,
  });

  $runner->init;
  # Force setting the package variables as the 'run' method is not used here
  if ($set_package_variables) {
    $runner->_set_package_variables();
  }

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
