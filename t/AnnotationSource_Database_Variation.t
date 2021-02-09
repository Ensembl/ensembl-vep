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

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation');


## DATABASE TESTS
#################

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 25 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens') if $can_use_db;

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Config');

  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
  });
  ok($cfg, 'get new config object');
  
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
    config => $cfg
  });

  ok($as, 'new is defined');

  is(ref($as), 'Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation', 'check class');

  throws_ok {
    Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
      config => Bio::EnsEMBL::VEP::Config->new({
        %$db_cfg,
        database => 1,
        offline => 0,
        species => 'homo_vepiens',
        check_frequency => 1,
      })
    })
  } qr/not supported using database/, 'check_frequency not supported';

  $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
    config => $cfg
  });

  ## METHOD TESTS
  ###############

  is_deeply(
    $as->info,
    {
      'COSMIC' => '67',
      'dbSNP' => '138',
      'ClinVar' => '23', 
    },
    'info'
  );

  is($as->have_pubmed, 0, 'have_pubmed');
  is($as->phenotype_attrib_id, '418', 'phenotype_attrib_id');
  is($as->clinvar_source_id_cache, '32', 'clinvar_source_id_cache');
  
  ## TESTS WITH AN INPUT BUFFER
  #############################

  use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
  my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
  ok($p, 'get parser object');

  use_ok('Bio::EnsEMBL::VEP::InputBuffer');
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'check class');

  is(ref($ib->next()), 'ARRAY', 'check buffer next');

  is_deeply(
    $as->get_all_regions_by_InputBuffer($ib),
    [
      [
        '21',
        511
      ],
      [
        '21',
        512
      ],
      [
        '21',
        513
      ],
      [
        '21',
        514
      ],
      [
        '21',
        515
      ],
      [
        '21',
        517
      ],
      [
        '21',
        518
      ],
      [
        '21',
        519
      ]
    ],
    'get_all_regions_by_InputBuffer'
  );

  my $features = $as->get_all_features_by_InputBuffer($ib);
  is(ref($features), 'ARRAY', 'get_all_features_by_InputBuffer ref 1');
  is(ref($features->[0]), 'HASH', 'get_all_features_by_InputBuffer ref 2');
  is(ref($features->[-1]), 'HASH', 'get_all_features_by_InputBuffer ref 3');
  is($features->[0]->{variation_name}, 'rs142513484', 'get_all_features_by_InputBuffer variation_name');
  is(scalar @$features, 133, 'get_all_features_by_InputBuffer count');

  # do it again to get them from memory
  $features = $as->get_all_features_by_InputBuffer($ib);
  is($features->[0]->{variation_name}, 'rs142513484', 'get_all_features_by_InputBuffer again');

  $ib->next();
  is_deeply($as->get_all_features_by_InputBuffer($ib), [], 'get_all_features_by_InputBuffer on empty buffer');

  # reset
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  $as->annotate_InputBuffer($ib);
  my $vf = $ib->buffer->[0];

  is_deeply(
    $vf->{existing},
    [{
      'phenotype_or_disease' => 0,
      'failed' => 0,
      'somatic' => 0,
      'strand' => 1,
      'allele_string' => 'C/T',
      'minor_allele_freq' => '0.000998403',
      'clin_sig' => undef,
      'end' => 25585733,
      'variation_name' => 'rs142513484',
      'minor_allele' => 'T',
      'start' => 25585733,
      'variation_id' => 28751744,
      'matched_alleles' => [{
        a_allele => 'T',
        a_index => 0,
        b_allele => 'T',
        b_index => 0,
      }]
    }],
    'annotate_InputBuffer'
  );

  is(scalar (grep {$_->{existing}} @{$ib->buffer}), 132, 'annotate_InputBuffer count annotated');

  $vf->{start}++;
  delete $vf->{existing};
  $ib->{buffer} = [$vf];
  $as->annotate_InputBuffer($ib);

  is($vf->{existing}, undef, 'annotate_InputBuffer - miss by one');


  # construct one to test phenotype_or_disease and clin_sig
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->create_input_file([qw(21 25891796 . C T . . .)])});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next;

  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{existing},
    [{
      'phenotype_or_disease' => 1,
      'failed' => 0,
      'somatic' => 0,
      'strand' => 1,
      'allele_string' => 'C/T',
      'minor_allele_freq' => '0.000199681',
      'clin_sig' => 'not_provided,pathogenic',
      'end' => 25891796,
      'variation_name' => 'rs63750066',
      'minor_allele' => 'T',
      'variation_id' => 13264416,
      'start' => 25891796,
      'matched_alleles' => [{
        a_allele => 'T',
        a_index => 0,
        b_allele => 'T',
        b_index => 0,
      }]
    }],
    'annotate_InputBuffer - phenotype_or_disease'
  );

  1;
}


# done
done_testing();
