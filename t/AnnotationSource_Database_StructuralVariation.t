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
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation');


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
  skip 'No local database configured', 19 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens') if $can_use_db;

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Config');

  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    everything => 1,
    xref_refseq => 1,
  });
  ok($cfg, 'get new config object');
  
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation->new({
    config => $cfg
  });

  ok($as, 'new is defined');

  is(ref($as), 'Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation', 'check class');


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
  is(ref($features->[0]), 'Bio::EnsEMBL::Variation::StructuralVariationFeature', 'get_all_features_by_InputBuffer ref 2');
  is(ref($features->[-1]), 'Bio::EnsEMBL::Variation::StructuralVariationFeature', 'get_all_features_by_InputBuffer ref 3');
  is($features->[0]->{variation_name}, 'nsv1061577', 'get_all_features_by_InputBuffer variation_name');
  is(scalar @$features, 750, 'get_all_features_by_InputBuffer count');

  # do it again to get them from memory
  $features = $as->get_all_features_by_InputBuffer($ib);
  is($features->[0]->{variation_name}, 'nsv1061577', 'get_all_features_by_InputBuffer again');

  $ib->next();
  is_deeply($as->get_all_features_by_InputBuffer($ib), [], 'get_all_features_by_InputBuffer on empty buffer');

  # reset
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  $as->annotate_InputBuffer($ib);
  my $vf = $ib->buffer->[0];

  is_deeply(
    $vf->{overlapping_svs},
    {
      'nsv1189327' => 1,
      'nsv531521' => 1,
      'nsv996057' => 1,
      'nsv916276' => 1,
      'nsv1061577' => 1,
      'nsv544385' => 1,
      'nsv1055407' => 1,
      'nsv544410' => 1,
      'nsv544386' => 1
    },
    'annotate_InputBuffer'
  );

  is(scalar (grep {$_->{overlapping_svs}} @{$ib->buffer}), 132, 'annotate_InputBuffer count annotated');

  1;
}


# done
done_testing();
