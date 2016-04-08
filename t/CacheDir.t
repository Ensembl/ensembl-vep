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
use_ok('Bio::EnsEMBL::VEP::CacheDir');

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg_hash = $test_cfg->base_testing_cfg;

my $cfg = Bio::EnsEMBL::VEP::Config->new($cfg_hash);
ok($cfg, 'get new config object');

my $cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}});
ok($cd, 'new is defined');

is(ref($cd), 'Bio::EnsEMBL::VEP::CacheDir', 'check class');

# should be able to work out the dir without species (uses default) or assembly (scans dir)
ok(
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({
      offline => 1,
      cache_version => $test_cfg->{cache_version}
    })
  }),
  "new with only version"
);

# and should be able to pass full path
ok(
  Bio::EnsEMBL::VEP::CacheDir->new({
    dir => $test_cfg->{cache_dir},
    config => Bio::EnsEMBL::VEP::Config->new({offline => 1})
  }),
  "new with full dir path"
);


# various fails on non-existent dir
throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg})
} qr/No root_dir or dir specified/, 'new without root_dir';

throws_ok { 
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, species => 'foo'})
  })
} qr/Cache directory .+ not found/, 'new with invalid species';

throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, cache_version => 20})
  })
} qr/No cache found for .+ 20/, 'new with invalid version';

throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, assembly => 'bar'})
  })
} qr/Cache assembly version .+bar.+ do not match/, 'new with invalid assembly';


## METHOD CALLS
###############

$cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}});

is($cd->root_dir, $cfg_hash->{dir}, 'root_dir');
is($cd->dir, $test_cfg->{cache_dir}, 'dir');

is_deeply(
  $cd->info,
  {
    'polyphen' => 'b',
    'sift' => 'b',
    'polyphen_version' => '2.2.2',
    'species' => 'homo_sapiens',
    'version_data' => {
      'polyphen' => '2.2.2',
      'sift' => 'sift5.2.2',
      'COSMIC' => '71',
      'ESP' => '20140509',
      'gencode' => 'GENCODE',
      'HGMD-PUBLIC' => '20142',
      'genebuild' => '2014-07',
      'regbuild' => '13.0',
      'assembly' => 'GRCh38.p2',
      'dbSNP' => '138',
      'ClinVar' => '201410'
    },
    'build' => 'all',
    'port' => '3306',
    'regulatory' => '1',
    'host' => 'genebuild13',
    'sift_version' => 'sift5.0.2',
    'cell_types' => 'HeLa-S3,GM06990,U2OS,CD4,IMR90,HL-60,HepG2,Lymphoblastoid,CD133,CD36,K562,GM12878,HUVEC,NHEK,H1ESC,MultiCell,K562b,NH-A,HSMM,HMEC,A549,AG04449,AG04450,AG09309,AG09319,AG10803,Caco-2,Chorion,CMK,GM10847,GM12801,GM12864,GM12865,GM12872,GM12873,GM12874,GM12875,GM12891,GM12892,GM15510,GM18505,GM18507,GM18526,GM18951,GM19099,GM19193,GM19238,GM19239,GM19240,H7ESC,H9ESC,HAEpiC,HCF,HCM,HCPEpiC,HCT116,HEEpiC,HEK293b,HEK293,HepG2b,HGF,HIPEpiC,HNPCEpiC,HRCEpiC,HRE,HRPEpiC,Jurkat,LHSR,MCF7,Medullo,Melano,NB4,NHBE,NHDF-neo,NHLF,NT2-D1,Panc1,PanIslets,PFSK1,SAEC,SKMC,SKNMC,SKNSHRA,Th1,Th2,WERIRB1,RPTEC,ProgFib,HSMMtube,Osteobl,MCF10A-Er-Src,HPAEpiC,Fibrobl,GM12878-XiMat,BJ,NHDF-AD,Monocytes-CD14+,DND-41',
    'variation_cols' => [
      'variation_name',
      'failed',
      'somatic',
      'start',
      'end',
      'allele_string',
      'strand',
      'minor_allele',
      'minor_allele_freq',
      'clin_sig',
      'phenotype_or_disease',
      'pubmed',
      'AFR',
      'AMR',
      'EAS',
      'EUR',
      'SAS',
      'AA',
      'EA',
      'ExAC',
      'ExAC_AFR',
      'ExAC_AMR',
      'ExAC_Adj',
      'ExAC_EAS',
      'ExAC_FIN',
      'ExAC_NFE',
      'ExAC_OTH',
      'ExAC_SAS'
    ],
    'user' => 'ensro',
    'assembly' => 'GRCh38'
  },
  'info'
);

is_deeply(
  $cd->version_data,
  {
    'polyphen' => '2.2.2',
    'sift' => 'sift5.2.2',
    'COSMIC' => '71',
    'ESP' => '20140509',
    'gencode' => 'GENCODE',
    'HGMD-PUBLIC' => '20142',
    'genebuild' => '2014-07',
    'regbuild' => '13.0',
    'assembly' => 'GRCh38.p2',
    'dbSNP' => '138',
    'ClinVar' => '201410'
  },
  'version_data'
);

my $as;
ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources');

is(ref($as), 'ARRAY', 'get_all_AnnotationSources return type');
is(scalar @$as, 1, 'get_all_AnnotationSources count');
is(ref($as->[0]), 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript', 'first AnnotationSource is Transcript type');

# switch on regulatory
$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, regulatory => 1});
ok($cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}}), 'new with regulatory');

ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources with regulatory');
is(scalar @$as, 2, 'get_all_AnnotationSources with regulatory count');
is(ref($as->[1]), 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat', 'second AnnotationSource is RegFeat type');

done_testing();
