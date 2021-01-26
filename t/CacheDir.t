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
use File::Copy;

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

# copy chr_synonyms.txt file temporarily
copy($test_cfg->{chr_synonyms}, $cd->dir);

# detects FASTA file and synonyms
$cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}});
is($cd->param('synonyms'), $cd->dir.'/chr_synonyms.txt', 'detect synonyms');
ok($cd->param('fasta') =~ /test\.fa/, 'detect FASTA');

# don't bump off existing entries
$cd = Bio::EnsEMBL::VEP::CacheDir->new({
  config => Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, fasta => 'foo', synonyms => 'bar'}),
  root_dir => $cfg_hash->{dir}
});
is($cd->param('synonyms'), 'bar', 'dont overwrite synonyms param');
is($cd->param('fasta'), 'foo', 'dont overwrite fasta param');

unlink($cd->dir.'/chr_synonyms.txt');

$cd = Bio::EnsEMBL::VEP::CacheDir->new({
  root_dir => $cfg_hash->{dir},
  config => Bio::EnsEMBL::VEP::Config->new({
    %$cfg_hash,
    no_fasta => 1,
  })
});
ok(!$cd->param('fasta'), 'switch off detect FASTA');


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
    config => Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, species => 'bar'})
  })->dir()
} qr/Cache directory .+ not found/, 'new with invalid species';

throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, cache_version => 20})
  })->dir()
} qr/No cache found for .+ 20/, 'new with invalid version';

throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, assembly => 'bar'})
  })->dir()
} qr/Cache assembly version .+bar.+ do not match/, 'new with invalid assembly';


my %tmp_cfg_hash = %$cfg_hash;
delete $tmp_cfg_hash{assembly};

throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%tmp_cfg_hash, species => 'foo'})
  })->dir()
} qr/Multiple assemblies found/, 'new with no assembly, multiple available';

throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%tmp_cfg_hash, species => 'foo', assembly => 'tar'})
  })->dir()
} qr/No cache found/, 'new with assembly, multiple available';

throws_ok {
  Bio::EnsEMBL::VEP::CacheDir->new({
    root_dir => $cfg_hash->{dir},
    config => Bio::EnsEMBL::VEP::Config->new({%tmp_cfg_hash, species => 'foo', assembly => 'bar'})
  })->info()
} qr/Mismatch in assembly versions/, 'new, mismatch assembly';


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
    'species' => 'homo_sapiens',
    'version_data' => {
      'sift' => 'sift5.2.2',
      'polyphen' => '2.2.2',
      '1000genomes' => 'phase3',
      'COSMIC' => '80',
      'ESP' => 'V2-SSA137',
      'gnomAD' => '170228',
      'gencode' => 'GENCODE 24',
      'genebuild' => '2014-07',
      'HGMD-PUBLIC' => '20164',
      'regbuild' => '13.0',
      'dbSNP' => '149',
      'ClinVar' => '201704',
      'assembly' => 'GRCh38.p5'
    },
    'build' => 'all',
    'regulatory' => '1',
    'cell_types' => 'A549,CD14+CD16-Monocyte:CordBlood:hist,CD14+CD16-Monocyte:VenousBlood:hist,CD4+abTcell:VenousBlood:hist,CD8+abTcell:CordBlood:hist,CMCD4+abTcell:VenousBlood:hist,DND-41,eosinophil:VenousBlood:hist,EPC:VenousBlood:hist,erythroblast:CordBlood:hist,GM12878,H1ESC,HeLa-S3,HepG2,HMEC,HSMM,HSMMtube,HUVEC,HUVECprol:CordBlood:hist,IMR90,K562,M0Macrophage:CordBlood:hist,M0Macrophage:VenousBlood:hist,M1Macrophage:CordBlood:hist,M1Macrophage:VenousBlood:hist,M2Macrophage:CordBlood:hist,M2Macrophage:VenousBlood:hist,Monocytes-CD14+,MSC:VenousBlood:hist,naiveBcell:VenousBlood:hist,neutroMyelocyte:BoneMarrow:hist,neutrophil:CordBlood:hist,neutrophil:VenousBlood:hist,NH-A,NHDF-AD,NHEK,NHLF,Osteobl,MultiCell',
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
      'gnomAD',
      'gnomAD_AFR',
      'gnomAD_AMR',
      'gnomAD_ASJ',
      'gnomAD_EAS',
      'gnomAD_FIN',
      'gnomAD_NFE',
      'gnomAD_OTH',
      'gnomAD_SAS',
      'var_synonyms'
    ],
    'assembly' => 'GRCh38',
    'valid_chromosomes' => [21, 22, 'LRG_485'],
  },
  'info'
);

is_deeply(
  $cd->version_data,
  {
    'sift' => 'sift5.2.2',
    'polyphen' => '2.2.2',
    '1000genomes' => 'phase3',
    'COSMIC' => '80',
    'ESP' => 'V2-SSA137',
    'gnomAD' => '170228',
    'gencode' => 'GENCODE 24',
    'genebuild' => '2014-07',
    'HGMD-PUBLIC' => '20164',
    'regbuild' => '13.0',
    'dbSNP' => '149',
    'ClinVar' => '201704',
    'assembly' => 'GRCh38.p5'
  },
  'version_data'
);

my $as;
ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources');

is(ref($as), 'ARRAY', 'get_all_AnnotationSources return type');
is(scalar @$as, 1, 'get_all_AnnotationSources count');
is(ref($as->[0]), 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript', 'first AnnotationSource is Transcript type');
is_deeply($as->[0]->info, $cd->version_data, 'version_data added to annotation source');

# switch on regulatory
$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, regulatory => 1});
ok($cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}}), 'new with regulatory');

ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources with regulatory');
is(scalar @$as, 2, 'get_all_AnnotationSources with regulatory count');
is(ref($as->[1]), 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat', 'second AnnotationSource is RegFeat type');
is_deeply(
  $as->[1]->{available_cell_types},
  [
    'A549',
    'CD14+CD16-Monocyte:CordBlood:hist',
    'CD14+CD16-Monocyte:VenousBlood:hist',
    'CD4+abTcell:VenousBlood:hist',
    'CD8+abTcell:CordBlood:hist',
    'CMCD4+abTcell:VenousBlood:hist',
    'DND-41',
    'eosinophil:VenousBlood:hist',
    'EPC:VenousBlood:hist',
    'erythroblast:CordBlood:hist',
    'GM12878',
    'H1ESC',
    'HeLa-S3',
    'HepG2',
    'HMEC',
    'HSMM',
    'HSMMtube',
    'HUVEC',
    'HUVECprol:CordBlood:hist',
    'IMR90',
    'K562',
    'M0Macrophage:CordBlood:hist',
    'M0Macrophage:VenousBlood:hist',
    'M1Macrophage:CordBlood:hist',
    'M1Macrophage:VenousBlood:hist',
    'M2Macrophage:CordBlood:hist',
    'M2Macrophage:VenousBlood:hist',
    'Monocytes-CD14+',
    'MSC:VenousBlood:hist',
    'naiveBcell:VenousBlood:hist',
    'neutroMyelocyte:BoneMarrow:hist',
    'neutrophil:CordBlood:hist',
    'neutrophil:VenousBlood:hist',
    'NH-A',
    'NHDF-AD',
    'NHEK',
    'NHLF',
    'Osteobl',
    'MultiCell'
  ],
  'cell types added from info'
);


# switch on variation
$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, check_existing => 1});
ok($cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}}), 'new with variation');

ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources with variation');
is(scalar @$as, 2, 'get_all_AnnotationSources with variation count');
is(ref($as->[1]), 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation', 'second AnnotationSource is Variation type');


# check gnomad/exac switchup
$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, check_existing => 1, af_gnomad => 1});
$cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}});
ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources - gnomad available');

$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, check_existing => 1, af_gnomad => 1});
$cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $test_cfg->{exac_root_dir}});
throws_ok {$cd->get_all_AnnotationSources} qr/gnomad.+not available.+exac/i, 'get_all_AnnotationSources - gnomad not available';

$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, check_existing => 1, af_exac => 1});
$cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $test_cfg->{exac_root_dir}});
ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources - exac available');

$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, check_existing => 1, af_exac => 1});
$cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}});
throws_ok {$cd->get_all_AnnotationSources} qr/exac.+not available.+gnomad/i, 'get_all_AnnotationSources - exac not available';


# switch on both
$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, check_existing => 1, regulatory => 1});
ok($cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}}), 'new with both');
ok($as = $cd->get_all_AnnotationSources, 'get_all_AnnotationSources with both');
is_deeply([map {(split('::', ref($_)))[-1]} @{$cd->get_all_AnnotationSources}], [qw(Transcript RegFeat Variation)], 'get_all_AnnotationSources with both - types');

# tabix type
$cfg = Bio::EnsEMBL::VEP::Config->new({%$cfg_hash, check_existing => 1});
$cd = Bio::EnsEMBL::VEP::CacheDir->new({config => $cfg, root_dir => $cfg_hash->{dir}});

$cd->info->{var_type} = 'tabix';
is(
  ref($cd->get_all_AnnotationSources->[1]),
  'Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix',
  'VariationTabix type'
);



done_testing();
