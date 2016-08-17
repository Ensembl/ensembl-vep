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

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSourceAdaptor');

my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new();
ok($asa, 'new is defined');

is(ref($asa), 'Bio::EnsEMBL::VEP::AnnotationSourceAdaptor', 'check class');

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg_hash = $test_cfg->base_testing_cfg;

my $cfg = Bio::EnsEMBL::VEP::Config->new($cfg_hash);
ok($cfg, 'get new config object');

ok($asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $cfg}), 'new with config');



## METHOD CALLS
###############

my $exp = [
  bless( {
    '_config' => $cfg,
    'serializer_type' => undef,
    'dir' => $test_cfg->{cache_dir},
    'cache_region_size' => $cfg->param('cache_region_size'),
    'gencode_basic' => undef,
    'all_refseq' => undef,
    'source_type' => 'ensembl',
    'sift' => undef,
    'polyphen' => undef,
    'everything' => undef,
    'info' => {
      'polyphen' => '2.2.2',
      'sift' => 'sift5.2.2',
      'COSMIC' => '75',
      'ESP' => '20141103',
      'gencode' => 'GENCODE 24',
      'HGMD-PUBLIC' => '20154',
      'genebuild' => '2014-07',
      'regbuild' => '13.0',
      'assembly' => 'GRCh38.p5',
      'dbSNP' => '146',
      'ClinVar' => '201601'
    },
    'valid_chromosomes' => [21, 'LRG_485'],
  }, 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript' )
];

is_deeply($asa->get_all_from_cache(), $exp, 'get_all_from_cache');

is_deeply($asa->get_all(), $exp, 'get_all');

$asa->param('check_existing', 1);
is_deeply(ref($asa->get_all()->[0]), 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation', 'get_all - var comes first');
$asa->param('check_existing', 0);

$asa->param('custom', [$test_cfg->{custom_vcf}]);
throws_ok {$asa->get_all_custom} qr/No format/, 'get_all_custom - no format';

$asa->param('custom', [$test_cfg->{custom_vcf}.',test,foo,exact']);
throws_ok {$asa->get_all_custom} qr/Unknown or unsupported format foo/, 'get_all_custom - invalid format';

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 4 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

  $asa->param('custom', [$test_cfg->{custom_vcf}.',test,vcf,exact']);
  is_deeply(
    $asa->get_all_custom(),
    [
      bless( {
        'short_name' => 'test',
        '_config' => $asa->config,
        'report_coords' => 0,
        'file' => $test_cfg->{custom_vcf},
        'type' => 'exact',
        'info' => {
          'custom_info' => {
            'short_name' => 'test',
            'report_coords' => undef,
            'file' => $test_cfg->{custom_vcf},
            'type' => 'exact'
          }
        }
      }, 'Bio::EnsEMBL::VEP::AnnotationSource::File::VCF' )
    ],
    'get_all_custom'
  );

  $asa->param('custom', [$test_cfg->{custom_vcf}.',test,vcf']);
  is_deeply(
    $asa->get_all_custom(),
    [
      bless( {
        'short_name' => 'test',
        '_config' => $asa->config,
        'report_coords' => 0,
        'file' => $test_cfg->{custom_vcf},
        'type' => 'overlap',
        'info' => {
          'custom_info' => {
            'short_name' => 'test',
            'report_coords' => undef,
            'file' => $test_cfg->{custom_vcf},
            'type' => 'overlap'
          }
        }
      }, 'Bio::EnsEMBL::VEP::AnnotationSource::File::VCF' )
    ],
    'get_all_custom - default overlap type'
  );

  $asa->param('custom', [$test_cfg->{custom_vcf}.',test,vcf,overlap,1']);
  is_deeply(
    $asa->get_all_custom(),
    [
      bless( {
        'short_name' => 'test',
        '_config' => $asa->config,
        'report_coords' => 1,
        'file' => $test_cfg->{custom_vcf},
        'type' => 'overlap',
        'info' => {
          'custom_info' => {
            'short_name' => 'test',
            'report_coords' => 1,
            'file' => $test_cfg->{custom_vcf},
            'type' => 'overlap',
          }
        }
      }, 'Bio::EnsEMBL::VEP::AnnotationSource::File::VCF' )
    ],
    'get_all_custom - report_coords'
  );

  $asa->param('custom', [$test_cfg->{custom_vcf}.',test,vcf,overlap,1,FOO,BAR']);
  is_deeply(
    $asa->get_all_custom(),
    [
      bless( {
        'short_name' => 'test',
        '_config' => $asa->config,
        'report_coords' => 1,
        'fields' => ['FOO', 'BAR'],
        'file' => $test_cfg->{custom_vcf},
        'type' => 'overlap',
        'info' => {
          'custom_info' => {
            'short_name' => 'test',
            'report_coords' => 1,
            'fields' => ['FOO', 'BAR'],
            'file' => $test_cfg->{custom_vcf},
            'type' => 'overlap',
          }
        }
      }, 'Bio::EnsEMBL::VEP::AnnotationSource::File::VCF' )
    ],
    'get_all_custom - fields'
  );
}

$asa->param('custom', []);

## DATABASE TESTS
#################

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;
  my $can_use_db = $db_cfg && scalar keys %$db_cfg;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 5 unless $can_use_db;

  my $multi;

  if($can_use_db) {
    eval q{
      use Bio::EnsEMBL::Test::TestUtils;
      use Bio::EnsEMBL::Test::MultiTestDB;
      1;
    };

    $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens');
  }

  $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({
    config => Bio::EnsEMBL::VEP::Config->new({
      %$db_cfg,
      database => 1,
      offline => 0,
      species => 'homo_vepiens',
    })
  });

  $exp = [
    bless( {
      'polyphen' => undef,
      'sift' => undef,
      'everything' => undef,
      'uniprot' => undef,
      '_config' => $asa->config,
      'xref_refseq' => undef,
      'polyphen_analysis' => 'humvar',
      'cache_region_size' => 50000,
      'protein' => undef,
      'domains' => undef,
      'gencode_basic' => undef,
      'source_type' => 'ensembl',
      'core_type' => 'core',
      'all_refseq' => undef,
      'assembly' => undef,
      '_species' => 'homo_vepiens',
      'info' => {
        'genebuild' => '2014-07',
        'gencode' => 'GENCODE 24',
        'assembly' => 'GRCh38.p5'
      }
    }, 'Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript' )
  ];

  my $asas = $asa->get_all_from_database();
  delete $asas->[0]->{_adaptors};

  is_deeply($asas, $exp, 'get_all_from_database');

  $asas = $asa->get_all();
  delete $asas->[0]->{_adaptors};

  is_deeply($asas, $exp, 'get_all - DB');

  $asa->param('check_existing', 1);
  is_deeply(
    $asa->get_all_from_database()->[1],
    bless( {
      '_config' => $asa->config,
      'cache_region_size' => 50000,
      'failed' => 0,
      'no_check_alleles' => undef,
    }, 'Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation' ),
    'get_all_from_database - variation'
  );
  $asa->param('check_existing', 0);

  $asa->param('regulatory', 1);
  is_deeply(
    $asa->get_all_from_database()->[1],
    bless( {
      '_config' => $asa->config,
      'cache_region_size' => 50000,
      'cell_type' => undef,
    }, 'Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat' ),
    'get_all_from_database - regfeat'
  );
  $asa->param('regulatory', 0);

  $asa->param('check_svs', 1);
  is_deeply(
    $asa->get_all_from_database()->[1],
    bless( {
      '_config' => $asa->config,
      'cache_region_size' => 50000,
    }, 'Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation' ),
    'get_all_from_database - SV'
  );
  $asa->param('check_svs', 0);
};



done_testing();
