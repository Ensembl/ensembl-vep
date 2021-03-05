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

use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::StructuralVariationFeature;

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

my $cfg_hash = $test_cfg->base_testing_cfg;

use_ok('Bio::EnsEMBL::VEP::OutputFactory');

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'JSON module not available', 19 unless $Bio::EnsEMBL::VEP::OutputFactory::CAN_USE_JSON;

  ## BASIC TESTS
  ##############

  # use test
  use_ok('JSON');
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

  $vf->{existing}->[0]->{pubmed} = "10,20,30";
  $vf->{existing}->[0]->{clin_sig} = "pathogenic,benign";
  my $ex = $vf->{existing}->[0];

  my $frequency_hash = {Allele => 'T'};
  my $super_of = Bio::EnsEMBL::VEP::OutputFactory->new({config => $cfg});
  $super_of->{af_1kg} = 1;
  $super_of->{af_esp} = 1;
  $super_of->{af_gnomad} = 1;
  $super_of->{af_exac} = 1;
  $super_of->add_colocated_frequency_data({}, $frequency_hash, $ex);
  is_deeply(
    $of->add_colocated_variant_info_JSON({}, [$frequency_hash], $ex),
    {
      'colocated_variants' => [
        {
          'frequencies' => {
            'T' => {
              'amr' => '0.0014',
              'gnomad_sas' => '0',
              'gnomad' => '0.0003478',
              'ea' => '0',
              'gnomad_oth' => '0',
              'gnomad_asj' => '0',
              'gnomad_nfe' => '1.886e-05',
              'aa' => '0.004998',
              'gnomad_afr' => '0.004643',
              'afr' => '0.003',
              'gnomad_amr' => '0.0003236',
              'gnomad_fin' => '0',
              'sas' => '0',
              'gnomad_eas' => '0',
              'eur' => '0',
              'eas' => '0'
            }
          },
          'id' => 'rs142513484',
          'minor_allele_freq' => '0.0010',
          'minor_allele' => 'T',
          'end' => 25585733,
          'start' => 25585733,
          'strand' => 1,
          'allele_string' => 'C/T',
	  'var_synonyms' => {
	    'ClinVar' => ['TEST']
	  },
          'pubmed' => [10, 20, 30],
          'clin_sig' => ["pathogenic", "benign"],
        }
      ]
    },
    'add_colocated_variant_info_JSON'
  );

  # Frequency test for a variant that has a multi-allelic co-located variant.
  # The input allele matches the co-located variant allele with no associated
  # 1kg frequency.
  # 1kg frequency is availabe on the other allele of the co-located variant.
  # The test checks that 1kg frequency for the other allele is not reported.
  no warnings 'qw';
  $ib = get_annotated_buffer({
    input_file => $test_cfg->create_input_file([qw(21 25891785 var_1 G A . . .)]),
    check_existing => 1
  });

  $vf = $ib->buffer->[0];
  $ex = $vf->{existing}->[0];

  $frequency_hash = {Allele => 'A'};
  $super_of = Bio::EnsEMBL::VEP::OutputFactory->new({config => $cfg});
  $super_of->{af_1kg} = 1;
  $super_of->{af_esp} = 1;
  $super_of->{af_gnomad} = 1;
  $super_of->{af_exac} = 1;
  $super_of->add_colocated_frequency_data({}, $frequency_hash, $ex);
  my $x = $of->add_colocated_variant_info_JSON({}, [$frequency_hash], $ex);

  is_deeply(
    $of->add_colocated_variant_info_JSON({}, [$frequency_hash], $ex),
    {
      'colocated_variants' => [
         {
           'strand' => 1,
           'id' => 'rs145564988',
           'allele_string' => 'G/A/T',
           'minor_allele_freq' => '0.0016',
           'frequencies' => {
              'A' => {
                'gnomad_afr' => '6.534e-05',
                'gnomad_sas' => '3.249e-05',
                'gnomad_fin' => '0',
                'gnomad_amr' => '0',
                'gnomad' => '7.313e-05',
                'gnomad_oth' => '0',
                'gnomad_asj' => '0',
                'gnomad_eas' => '0',
                'gnomad_nfe' => '0.0001433'
              }
           },
           'end' => '25891785',
           'minor_allele' => 'T',
           'start' => '25891785'
         }
      ],
    },
    'add_colocated_variant_info_JSON_multi'
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
    $json->decode($lines[28])->{'regulatory_feature_consequences'},
    [{
      'consequence_terms' => [
        'regulatory_region_variant'
      ],
      'variant_allele' => 'G',
      'regulatory_feature_id' => 'ENSR00000140751',
      'impact' => 'MODIFIER',
      'biotype' => 'promoter'
    }],
    'get_all_lines_by_InputBuffer - regulatory information'
  );

  is_deeply(
    $json->decode($lines[0]),
    {
      'input' => "21\t25585733\trs142513484\tC\tT\t.\t.\t.\tGT\t0|0",
      'colocated_variants' => [
        {
          'frequencies' => {
            'T' => {
              'amr' => 0.0014,
              'gnomad_sas' => 0,
              'gnomad' => 0.0003478,
              'ea' => 0,
              'gnomad_oth' => 0,
              'gnomad_asj' => 0,
              'gnomad_nfe' =>  1.886e-05,
              'aa' => 0.004998,
              'gnomad_afr' => 0.004643,
              'afr' => 0.003,
              'gnomad_amr' => 0.0003236,
              'gnomad_fin' => 0,
              'sas' => 0,
              'gnomad_eas' => 0,
              'eur' => 0,
              'eas' => 0
            }
          },
          'id' => 'rs142513484',
          'minor_allele_freq' => 0.0010,
          'end' => 25585733,
          'minor_allele' => 'T',
          'start' => 25585733,
          'strand' => 1,
          'allele_string' => 'C/T',
          'var_synonyms' => {
            'ClinVar' => ['TEST']
          },
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
          'swissprot' => [
            'Q9NYK5'
          ],
          'hgvsc' => 'ENST00000307301.11:c.*18G>A',
          'hgnc_id' => 'HGNC:14027',
          'strand' => -1,
          'gene_symbol' => 'MRPL39',
          'exon' => '11/11',
          'cdna_start' => 1122,
          'transcript_id' => 'ENST00000307301',
          'gene_id' => 'ENSG00000154719',
          'canonical' => 1,
          'protein_id' => 'ENSP00000305682',
          'uniparc' => [
            'UPI00001AEAC0'
          ],
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
          'polyphen_score' => '0.001',
          'codons' => 'Gca/Aca',
          'swissprot' => [
            'Q9NYK5'
          ],
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
          'appris' => 'P1',
          'sift_prediction' => 'tolerated_low_confidence',
          'protein_id' => 'ENSP00000284967',
          'polyphen_prediction' => 'benign',
          'protein_start' => 331,
          'uniparc' => [
            'UPI00001AEE66'
          ],
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
          'gene_symbol_source' => 'Clone_based_ensembl_gene',
          'consequence_terms' => [
            'upstream_gene_variant'
          ],
          'strand' => -1,
          'gene_symbol' => 'AP000223.1',
          'transcript_id' => 'ENST00000567517',
          'impact' => 'MODIFIER'
        }
      ],
      'strand' => 1,
      'id' => 'rs142513484',
      'allele_string' => 'C/T',
      'most_severe_consequence' => 'missense_variant',
      'start' => 25585733,
    },
    'get_all_lines_by_InputBuffer - everything'
  );

  # check rejoin on minimal
  no warnings 'qw';
  $ib = get_annotated_buffer({
    input_file => $test_cfg->create_input_file([qw(21 25741665 . CAGAAGAAAG TAGAAGAAAG,C . . .)]),
    minimal => 1,
    pick_allele => 1,
  });
  $of = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $ib->config});

  is(scalar @{$ib->buffer}, 2, 'minimal - expanded count');
  is($ib->buffer->[0]->allele_string, 'C/T', 'minimal - expanded first allele string');

  $of->rejoin_variants_in_InputBuffer($ib);

  is(scalar @{$ib->buffer}, 1, 'minimal - rejoined count');
  is($ib->buffer->[0]->allele_string, 'CAGAAGAAAG/TAGAAGAAAG/C', 'minimal - rejoined allele string');

  is_deeply(
    {map {$_->{variant_allele} => $_} @{$json->decode($of->get_all_lines_by_InputBuffer($ib)->[0])->{transcript_consequences}}},
    {
      '-' => {
        'cds_start' => 67,
        'gene_id' => 'ENSG00000154727',
        'variant_allele' => '-',
        'cdna_end' => 603,
        'protein_start' => 23,
        'codons' => 'aAGAAGAAAGgc/agc',
        'cds_end' => 76,
        'consequence_terms' => [
          'inframe_deletion',
          'splice_region_variant'
        ],
        'protein_end' => 26,
        'strand' => 1,
        'amino_acids' => 'KKKG/S',
        'cdna_start' => 594,
        'transcript_id' => 'ENST00000354828',
        'impact' => 'MODERATE',
        'allele_num' => 2,
      },
      'T' => {
        'cds_start' => 67,
        'gene_id' => 'ENSG00000154727',
        'variant_allele' => 'T',
        'cdna_end' => 603,
        'protein_start' => 23,
        'codons' => 'Cca/Tca',
        'cds_end' => 76,
        'consequence_terms' => [
          'missense_variant'
        ],
        'protein_end' => 26,
        'strand' => 1,
        'amino_acids' => 'P/S',
        'cdna_start' => 594,
        'transcript_id' => 'ENST00000354828',
        'impact' => 'MODERATE',
        'allele_num' => 1,
      }
    },
    'minimal - get_all_lines_by_InputBuffer'
  );



  # test custom
  use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');

  SKIP: {
    no warnings 'once';

    ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
    skip 'Bio::DB::HTS::Tabix module not available', 2 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

    $ib = get_annotated_buffer({
      input_file => $test_cfg->{test_vcf},
      everything => 1,
      dir => $test_cfg->{cache_root_dir},
      custom => [$test_cfg->{custom_vcf}.',test,vcf,exact,,FOO'],
    });
    $of = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $ib->config});
    @lines = @{$of->get_all_lines_by_InputBuffer($ib)};

    is_deeply(
      $json->decode($lines[0])->{custom_annotations},
      {
        "test" => [
          {
            "fields" => {
              "FOO" => "BAR",
              "FILTER" => "PASS"
            },
            "name" => "test1",
            "allele" => "T"
          }
        ]
      },
      'custom_annotations'
    );

    $ib = get_annotated_buffer({
      input_data => "21\t25585733\t.\tCATG\tTACG",
      everything => 1,
      dir => $test_cfg->{cache_root_dir},
      custom => [$test_cfg->{custom_vcf}.',test,vcf,overlap'],
    });
    $of = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $ib->config});
    @lines = @{$of->get_all_lines_by_InputBuffer($ib)};

    is_deeply(
      $json->decode($lines[0])->{custom_annotations},
      {
        "test" => [
          {
            "name" => "test1",
          },
          {
            "name" => "del1",
          },
          {
            "name" => "del2",
          }

        ]
      },
      'custom_annotations overlap'
    );
  }
}

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
