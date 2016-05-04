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
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation');

my $dir = $test_cfg->{cache_dir};

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
ok($cfg, 'get new config object');

my $c = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation->new({
  config => $cfg,
  dir => $dir,
  cache_region_size => 1000000,
  cols => $test_cfg->{var_cols},
});
ok($c, 'new is defined');


## METHODS
##########

is($c->get_dump_file_name(1, '1-100'), $dir.'/1/1-100_var.gz', 'get_dump_file_name');
is($c->get_dump_file_name(1, 1, 100), $dir.'/1/1-100_var.gz', 'get_dump_file_name with end');

throws_ok { $c->get_dump_file_name() } qr/No chromosome/, 'get_dump_file_name no chromosome';
throws_ok { $c->get_dump_file_name(1) } qr/No region/, 'get_dump_file_name no region';

is($c->delimiter, qr/ /, 'delimiter');


# deserialization
my $features = $c->read_variations_from_file(
  $c->get_dump_file_name(
    $test_cfg->{cache_chr},
    $test_cfg->{cache_region}
  )
);
is(ref($features), 'ARRAY', 'read_variations_from_file ref');
is_deeply($features->[0], {
  'phenotype_or_disease' => 0,
  'somatic' => 0,
  'strand' => 1,
  'variation_name' => 'rs753123870',
  'allele_string' => 'G/C',
  'failed' => 0,
  'end' => '25973491',
  'start' => '25973491'
}, 'read_variations_from_file first el');

$features = $c->get_features_by_regions_uncached([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);
is(ref($features), 'ARRAY', 'get_features_by_regions_uncached ref 1');
is(ref($features->[0]), 'HASH', 'get_features_by_regions_uncached ref 2');
is($features->[0]->{variation_name}, 'rs753123870', 'get_features_by_regions_uncached variation_name');

# # now we should be able to retrieve the same from memory
$features = $c->get_features_by_regions_cached([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);
is(ref($features), 'ARRAY', 'get_features_by_regions_cached ref 1');
is(ref($features->[0]), 'HASH', 'get_features_by_regions_cached ref 2');
is($features->[0]->{variation_name}, 'rs753123870', 'get_features_by_regions_cached variation_name');

$c->clean_cache();
is_deeply($c->cache, {}, 'clean_cache');



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
  $c->get_all_regions_by_InputBuffer($ib),
  [[21, 25]],
  'get_all_regions_by_InputBuffer'
);

$features = $c->get_all_features_by_InputBuffer($ib);
is(ref($features), 'ARRAY', 'get_all_features_by_InputBuffer ref 1');
is(ref($features->[0]), 'HASH', 'get_all_features_by_InputBuffer ref 2');
is($features->[0]->{variation_name}, 'rs753123870', 'get_all_features_by_InputBuffer variation_name 1');
is($features->[-1]->{variation_name}, 'rs773346245', 'get_all_features_by_InputBuffer variation_name 2');
is(scalar @$features, 23650, 'get_all_features_by_InputBuffer count');

# do it again to get them from memory
$features = $c->get_all_features_by_InputBuffer($ib);
is($features->[0]->{variation_name}, 'rs753123870', 'get_all_features_by_InputBuffer again');

$ib->next();
is_deeply($c->get_all_features_by_InputBuffer($ib), [], 'get_all_features_by_InputBuffer on empty buffer');

# reset
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
$ib->next();

$c->annotate_InputBuffer($ib);
my $vf = $ib->buffer->[0];

is_deeply($vf->{existing}, [
  {
    'phenotype_or_disease' => 0,
    'ExAC_AMR' => 'T:0.000173',
    'SAS' => 'T:0.0000',
    'failed' => 0,
    'ExAC_NFE' => 'T:0',
    'AA' => 'T:0.005',
    'somatic' => 0,
    'ExAC_SAS' => 'T:0',
    'AFR' => 'T:0.0030',
    'strand' => 1,
    'allele_string' => 'C/T',
    'ExAC_Adj' => 'T:0.0004133',
    'minor_allele_freq' => '0.0010',
    'ExAC_FIN' => 'T:0',
    'AMR' => 'T:0.0014',
    'EUR' => 'T:0.0000',
    'clin_sig' => '',
    'EAS' => 'T:0.0000',
    'end' => 25585733,
    'ExAC' => 'T:4.119e-04',
    'ExAC_OTH' => 'T:0',
    'ExAC_AFR' => 'T:0.004681',
    'variation_name' => 'rs142513484',
    'ExAC_EAS' => 'T:0',
    'minor_allele' => 'T',
    'EA' => 'T:0',
    'start' => 25585733,
    'pubmed' => ''
  }
], 'annotate_InputBuffer');

is(scalar (grep {$_->{existing}} @{$ib->buffer}), 132, 'annotate_InputBuffer count annotated');

$vf->{start}++;
delete $vf->{existing};
$ib->{buffer} = [$vf];
$c->annotate_InputBuffer($ib);

is($vf->{existing}, undef, 'annotate_InputBuffer - miss by one');

# construct one to test phenotype_or_disease and clin_sig
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->create_input_file([qw(21 25891796 . C T . . .)])});
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
$ib->next;

$c->annotate_InputBuffer($ib);

is_deeply(
  $ib->buffer->[0]->{existing},
  [
    {
      'phenotype_or_disease' => '1',
      'ExAC_AMR' => 'T:0.0003457',
      'SAS' => 'T:0.0000',
      'failed' => 0,
      'ExAC_NFE' => 'T:2.998e-05',
      'AA' => undef,
      'somatic' => 0,
      'ExAC_SAS' => 'T:0',
      'AFR' => 'T:0.0000',
      'strand' => 1,
      'allele_string' => 'C/T',
      'ExAC_Adj' => 'T:5.768e-05',
      'minor_allele_freq' => '0.0002',
      'ExAC_FIN' => 'T:0',
      'AMR' => 'T:0.0014',
      'EUR' => 'T:0.0000',
      'clin_sig' => 'not_provided,pathogenic',
      'EAS' => 'T:0.0000',
      'end' => 25891796,
      'ExAC' => 'T:5.765e-05',
      'ExAC_OTH' => 'T:0.001101',
      'ExAC_AFR' => 'T:0',
      'variation_name' => 'rs63750066',
      'ExAC_EAS' => 'T:0',
      'minor_allele' => 'T',
      'EA' => undef,
      'start' => 25891796,
      'pubmed' => ''
    },
    {
      'phenotype_or_disease' => '1',
      'minor_allele_freq' => '',
      'clin_sig' => '',
      'failed' => 0,
      'end' => 25891796,
      'somatic' => 0,
      'strand' => 1,
      'variation_name' => 'CM930033',
      'allele_string' => 'HGMD_MUTATION',
      'minor_allele' => '',
      'start' => 25891796
    }
  ],
  'annotate_InputBuffer - phenotype_or_disease'
);

$p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
$ib->next();

$c->{old_maf} = 1;
$c->clean_cache();
$c->annotate_InputBuffer($ib);
is($ib->buffer->[0]->{existing}->[0]->{AMR}, 0.0014, 'old_maf');
$c->{old_maf} = 0;

# done
done_testing();
