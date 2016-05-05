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

use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::StructuralVariationFeature;

use lib $Bin;
use JSON;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

my $cfg_hash = $test_cfg->base_testing_cfg;

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::OutputFactory::JSON');

use_ok('Bio::EnsEMBL::VEP::Config');
use_ok('Bio::EnsEMBL::VEP::Runner');
my $cfg = Bio::EnsEMBL::VEP::Config->new();

my $of = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $cfg});

is(ref($of), 'Bio::EnsEMBL::VEP::OutputFactory::JSON', 'check class');



## METHOD TESTS
###############

is_deeply($of->headers, [], 'headers - empty');

my $ib = get_annotated_buffer({input_file => $test_cfg->{test_vcf}, check_existing => 1});

my $vf = $ib->buffer->[0];

is_deeply(
  $of->add_VariationFeatureOverlapAllele_info($vf, {}),
  {
    'transcript_consequences' => [
      {
        'gene_id' => 'ENSG00000154719',
        'variant_allele' => 'T',
        'cdna_end' => '1122',
        'consequence_terms' => [
          '3_prime_UTR_variant'
        ],
        'strand' => -1,
        'transcript_id' => 'ENST00000307301',
        'cdna_start' => '1122',
        'impact' => 'MODIFIER'
      },
      {
        'gene_id' => 'ENSG00000154719',
        'cds_start' => '991',
        'variant_allele' => 'T',
        'cdna_end' => '1033',
        'protein_start' => '331',
        'codons' => 'Gca/Aca',
        'cds_end' => '991',
        'consequence_terms' => [
          'missense_variant'
        ],
        'protein_end' => '331',
        'amino_acids' => 'A/T',
        'strand' => -1,
        'transcript_id' => 'ENST00000352957',
        'cdna_start' => '1033',
        'impact' => 'MODERATE'
      },
      {
        'gene_id' => 'ENSG00000260583',
        'variant_allele' => 'T',
        'distance' => 2407,
        'consequence_terms' => [
          'upstream_gene_variant'
        ],
        'strand' => -1,
        'transcript_id' => 'ENST00000567517',
        'impact' => 'MODIFIER'
      }
    ],
    'most_severe_consequence' => 'missense_variant'
  },
  'add_VariationFeatureOverlapAllele_info'
);

is_deeply(
  $of->add_colocated_variant_info($vf, {}),
  {
    'colocated_variants' => [
      {
        'aa_maf' => '0.005',
        'exac_adj_allele' => 'T',
        'exac_sas_maf' => '0',
        'ea_maf' => '0',
        'exac_fin_allele' => 'T',
        'exac_eas_allele' => 'T',
        'exac_afr_maf' => '0.004681',
        'eas_allele' => 'T',
        'amr_maf' => '0.0014',
        'exac_oth_allele' => 'T',
        'id' => 'rs142513484',
        'exac_eas_maf' => '0',
        'sas_allele' => 'T',
        'sas_maf' => '0.0000',
        'amr_allele' => 'T',
        'exac_allele' => 'T',
        'minor_allele_freq' => '0.0010',
        'eas_maf' => '0.0000',
        'end' => 25585733,
        'exac_fin_maf' => '0',
        'exac_maf' => '4.119e-04',
        'exac_nfe_allele' => 'T',
        'eur_allele' => 'T',
        'ea_allele' => 'T',
        'minor_allele' => 'T',
        'exac_afr_allele' => 'T',
        'start' => 25585733,
        'exac_oth_maf' => '0',
        'exac_adj_maf' => '0.0004133',
        'exac_sas_allele' => 'T',
        'strand' => 1,
        'aa_allele' => 'T',
        'exac_amr_allele' => 'T',
        'allele_string' => 'C/T',
        'exac_nfe_maf' => '0',
        'afr_allele' => 'T',
        'exac_amr_maf' => '0.000173',
        'afr_maf' => '0.0030',
        'eur_maf' => '0.0000'
      }
    ]
  },
  'add_colocated_variant_info'
);

$ib = get_annotated_buffer({input_file => $test_cfg->{test_vcf}});
$of = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $ib->config});

my @lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is(scalar @lines, scalar @{$ib->buffer}, 'get_all_lines_by_InputBuffer - count');

# we need to decode the JSON otherwise we run into hash order issues
my $json = JSON->new; 

is_deeply(
  $json->decode($lines[0]),
  {
    'input' => "21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0",
    'assembly_name' => 'GRCh38',
    'end' => 25585733,
    'seq_region_name' => '21',
    'strand' => 1,
    'transcript_consequences' => [
      {
        'gene_id' => 'ENSG00000154719',
        'consequence_terms' => [
          '3_prime_UTR_variant'
        ],
        'variant_allele' => 'T',
        'strand' => -1,
        'cdna_end' => 1122,
        'cdna_start' => 1122,
        'transcript_id' => 'ENST00000307301',
        'impact' => 'MODIFIER'
      },
      {
        'cds_start' => 991,
        'gene_id' => 'ENSG00000154719',
        'variant_allele' => 'T',
        'cdna_end' => 1033,
        'protein_start' => 331,
        'codons' => 'Gca/Aca',
        'cds_end' => 991,
        'consequence_terms' => [
          'missense_variant'
        ],
        'protein_end' => 331,
        'strand' => -1,
        'amino_acids' => 'A/T',
        'cdna_start' => 1033,
        'transcript_id' => 'ENST00000352957',
        'impact' => 'MODERATE'
      },
      {
        'gene_id' => 'ENSG00000260583',
        'consequence_terms' => [
          'upstream_gene_variant'
        ],
        'distance' => 2407,
        'variant_allele' => 'T',
        'strand' => -1,
        'transcript_id' => 'ENST00000567517',
        'impact' => 'MODIFIER'
      }
    ],
    'id' => 'rs142513484',
    'allele_string' => 'C/T',
    'most_severe_consequence' => 'missense_variant',
    'start' => 25585733
  },
  'get_all_lines_by_InputBuffer - check first'
);


$ib = get_annotated_buffer({
  input_file => $test_cfg->{test_vcf},
  everything => 1,
  dir => $test_cfg->{cache_root_dir},
});
$of = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $ib->config});
@lines = @{$of->get_all_lines_by_InputBuffer($ib)};

is_deeply(
  $json->decode($lines[0]),
  {
    'input' => "21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0",
    'colocated_variants' => [
      {
        'aa_maf' => 0.005,
        'exac_adj_allele' => 'T',
        'ea_maf' => 0,
        'exac_sas_maf' => 0,
        'exac_fin_allele' => 'T',
        'exac_eas_allele' => 'T',
        'exac_afr_maf' => 0.004681,
        'eas_allele' => 'T',
        'amr_maf' => 0.0014,
        'exac_oth_allele' => 'T',
        'sas_allele' => 'T',
        'exac_eas_maf' => 0,
        'id' => 'rs142513484',
        'sas_maf' => 0,
        'amr_allele' => 'T',
        'exac_allele' => 'T',
        'minor_allele_freq' => 0.001,
        'eas_maf' => 0,
        'end' => 25585733,
        'exac_fin_maf' => 0,
        'exac_maf' => 0.0004119,
        'eur_allele' => 'T',
        'exac_nfe_allele' => 'T',
        'minor_allele' => 'T',
        'ea_allele' => 'T',
        'start' => 25585733,
        'exac_afr_allele' => 'T',
        'exac_oth_maf' => 0,
        'exac_adj_maf' => 0.0004133,
        'exac_sas_allele' => 'T',
        'strand' => 1,
        'aa_allele' => 'T',
        'allele_string' => 'C/T',
        'exac_amr_allele' => 'T',
        'afr_allele' => 'T',
        'exac_nfe_maf' => 0,
        'exac_amr_maf' => 0.000173,
        'afr_maf' => 0.003,
        'eur_maf' => 0
      },
      {
        'aa_maf' => 0.005,
        'exac_adj_allele' => 'T',
        'ea_maf' => 0,
        'exac_sas_maf' => 0,
        'exac_fin_allele' => 'T',
        'exac_eas_allele' => 'T',
        'exac_afr_maf' => 0.004681,
        'eas_allele' => 'T',
        'amr_maf' => 0.0014,
        'exac_oth_allele' => 'T',
        'sas_allele' => 'T',
        'exac_eas_maf' => 0,
        'id' => 'rs142513484',
        'sas_maf' => 0,
        'amr_allele' => 'T',
        'exac_allele' => 'T',
        'minor_allele_freq' => 0.001,
        'eas_maf' => 0,
        'end' => 25585733,
        'exac_fin_maf' => 0,
        'exac_maf' => 0.0004119,
        'eur_allele' => 'T',
        'exac_nfe_allele' => 'T',
        'minor_allele' => 'T',
        'ea_allele' => 'T',
        'start' => 25585733,
        'exac_afr_allele' => 'T',
        'exac_oth_maf' => 0,
        'exac_adj_maf' => 0.0004133,
        'exac_sas_allele' => 'T',
        'strand' => 1,
        'aa_allele' => 'T',
        'allele_string' => 'C/T',
        'exac_amr_allele' => 'T',
        'afr_allele' => 'T',
        'exac_nfe_maf' => 0,
        'exac_amr_maf' => 0.000173,
        'afr_maf' => 0.003,
        'eur_maf' => 0
      }
    ],
    'assembly_name' => 'GRCh38',
    'end' => 25585733,
    'seq_region_name' => '21',
    'variant_class' => 'SNV',
    'transcript_consequences' => [
      {
        'variant_allele' => 'T',
        'cdna_end' => 1122,
        'swissprot' => 'Q9NYK5',
        'hgvsc' => 'ENST00000307301.11:c.*18G>A',
        'hgnc_id' => 'HGNC:14027',
        'strand' => -1,
        'gene_symbol' => 'MRPL39',
        'exon' => '11/11',
        'cdna_start' => 1122,
        'transcript_id' => 'ENST00000307301',
        'gene_id' => 'ENSG00000154719',
        'canonical' => 1,
        'appris' => 'A2',
        'protein_id' => 'ENSP00000305682',
        'uniparc' => 'UPI00001AEAC0',
        'biotype' => 'protein_coding',
        'gene_symbol_source' => 'HGNC',
        'consequence_terms' => [
          '3_prime_UTR_variant'
        ],
        'ccds' => 'CCDS33522.1',
        'impact' => 'MODIFIER',
        'tsl' => 5
      },
      {
        'hgvsp' => 'ENSP00000284967.6:p.Ala331Thr',
        'variant_allele' => 'T',
        'cdna_end' => 1033,
        'polyphen_score' => '0.021',
        'codons' => 'Gca/Aca',
        'swissprot' => 'Q9NYK5',
        'hgvsc' => 'ENST00000352957.8:c.991G>A',
        'protein_end' => 331,
        'strand' => -1,
        'amino_acids' => 'A/T',
        'hgnc_id' => 'HGNC:14027',
        'gene_symbol' => 'MRPL39',
        'exon' => '10/10',
        'transcript_id' => 'ENST00000352957',
        'cdna_start' => 1033,
        'gene_id' => 'ENSG00000154719',
        'cds_start' => 991,
        'appris' => 'P3',
        'sift_prediction' => 'tolerated_low_confidence',
        'protein_id' => 'ENSP00000284967',
        'polyphen_prediction' => 'benign',
        'protein_start' => 331,
        'uniparc' => 'UPI00001AEE66',
        'biotype' => 'protein_coding',
        'gene_symbol_source' => 'HGNC',
        'sift_score' => '0.17',
        'cds_end' => 991,
        'consequence_terms' => [
          'missense_variant'
        ],
        'ccds' => 'CCDS13573.1',
        'tsl' => 1,
        'impact' => 'MODERATE'
        },
      {
        'gene_id' => 'ENSG00000260583',
        'canonical' => 1,
        'distance' => 2407,
        'variant_allele' => 'T',
        'biotype' => 'antisense',
        'gene_symbol_source' => 'Clone_based_vega_gene',
        'consequence_terms' => [
          'upstream_gene_variant'
        ],
        'strand' => -1,
        'gene_symbol' => 'AP000223.42',
        'transcript_id' => 'ENST00000567517',
        'impact' => 'MODIFIER'
      }
    ],
    'strand' => 1,
    'id' => 'rs142513484',
    'regulatory_feature_consequences' => [
      {
        'consequence_terms' => [
          'regulatory_region_variant'
        ],
        'variant_allele' => 'T',
        'regulatory_feature_id' => 'ENSR00001963192',
        'impact' => 'MODIFIER',
        'biotype' => 'TF_binding_site'
      }
    ],
    'allele_string' => 'C/T',
    'most_severe_consequence' => 'missense_variant',
    'start' => 25585733
  },
  'get_all_lines_by_InputBuffer - everything'
);

# done
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
