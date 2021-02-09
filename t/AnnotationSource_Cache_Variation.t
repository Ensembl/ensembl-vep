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

ok($c->filter_variation({failed => 0}),  'filter_variation pass');
ok(!$c->filter_variation({failed => 1}), 'filter_variation fail');

$c->{failed} = 1;
ok($c->filter_variation({failed => 1}), 'filter_variation pass with failed on');
$c->{failed} = 0;


## NOVEL TESTS
##############

use_ok('Bio::EnsEMBL::Variation::VariationFeature');

my $input = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
  start  => 10,
  end    => 10,
  strand => 1,
  allele_string => 'A/G',
});

my $existing = {
  start  => 10,
  end    => 10,
  strand => 1,
  allele_string => 'A/G',
};

is_deeply(
  $c->compare_existing($input, $existing),
  {
    %$existing,
    matched_alleles => [{
      a_index => 0,
      a_allele => 'G',
      b_index => 0,
      b_allele => 'G',
    }]
  },
  'compare_existing exact match'
);

$input->{allele_string} = 'A/G/T';
is_deeply(
  $c->compare_existing($input, $existing),
  {
    %$existing,
    matched_alleles => [{
      a_index => 0,
      a_allele => 'G',
      b_index => 0,
      b_allele => 'G',
    }]
  },
  'compare_existing exact match from multiple'
);
$input->{allele_string} = 'A/G';

$existing->{allele_string} = 'A/T';
is(
  $c->compare_existing($input, $existing),
  undef,
  'compare_existing alleles dont match with check'
);

$c->{no_check_alleles} = 1;
is(
  $c->compare_existing($input, $existing),
  $existing,
  'compare_existing alleles dont match with check'
);
$c->{no_check_alleles} = 0;

$existing = {
  start  => 10,
  end    => 10,
  strand => -1,
  allele_string => 'T/C',
};

is_deeply(
  $c->compare_existing($input, $existing),
  {
    %$existing,
    matched_alleles => [{
      a_index => 0,
      a_allele => 'G',
      b_index => 0,
      b_allele => 'C',
    }]
  },
  'compare_existing rev strand exact match'
);

$existing->{allele_string} = 'T/G';
is(
  $c->compare_existing($input, $existing),
  undef,
  'compare_existing rev strand alleles dont match with check'
);

$c->{no_check_alleles} = 1;
is_deeply(
  $c->compare_existing($input, $existing),
  $existing,
  'compare_existing rev strand alleles dont match no check'
);
$c->{no_check_alleles} = 0;

# known variant missing alleles
$existing->{allele_string} = 'NULL';

is_deeply(
  $c->compare_existing($input, $existing),
  $existing,
  'compare_existing missing alleles'
);

$c->{no_check_alleles} = 1;
is_deeply(
  $c->compare_existing($input, $existing),
  $existing,
  'compare_existing missing alleles no check'
);
$c->{no_check_alleles} = 0;

$c->{exclude_null_alleles} = 1;
is($c->compare_existing($input, $existing), undef, 'compare_existing missing alleles exclude_null_alleles');
$c->{exclude_null_alleles} = 0;



## OTHER METHODS
################

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
  'SAS' => 'T:0',
  'failed' => 0,
  'somatic' => 0,
  'AFR' => 'T:0',
  'strand' => 1,
  'allele_string' => 'C/T',
  'minor_allele_freq' => '0.0002',
  'AMR' => 'T:0',
  'EUR' => 'T:0.001',
  'clin_sig' => '',
  'EAS' => 'T:0',
  'end' => '25000001',
  'variation_name' => 'rs574523538',
  'minor_allele' => 'T',
  'start' => '25000001',
  'pubmed' => ''
}, 'read_variations_from_file first el');

$features = $c->get_features_by_regions_uncached([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);
is(ref($features), 'ARRAY', 'get_features_by_regions_uncached ref 1');
is(ref($features->[0]), 'HASH', 'get_features_by_regions_uncached ref 2');
is($features->[0]->{variation_name}, 'rs574523538', 'get_features_by_regions_uncached variation_name');

# # now we should be able to retrieve the same from memory
$features = $c->get_features_by_regions_cached([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);
is(ref($features), 'ARRAY', 'get_features_by_regions_cached ref 1');
is(ref($features->[0]), 'HASH', 'get_features_by_regions_cached ref 2');
is($features->[0]->{variation_name}, 'rs574523538', 'get_features_by_regions_cached variation_name');

$c->clean_cache();
is_deeply($c->cache, {}, 'clean_cache');



## TESTS WITH AN INPUT BUFFER
#############################

use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
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
is($features->[0]->{variation_name}, 'rs142513484', 'get_all_features_by_InputBuffer variation_name 1');
is($features->[-1]->{variation_name}, 'COSM5057537', 'get_all_features_by_InputBuffer variation_name 2');
is(scalar @$features, 50848, 'get_all_features_by_InputBuffer count');

# do it again to get them from memory
$features = $c->get_all_features_by_InputBuffer($ib);
is($features->[0]->{variation_name}, 'rs142513484', 'get_all_features_by_InputBuffer again');

$ib->next();
is_deeply($c->get_all_features_by_InputBuffer($ib), [], 'get_all_features_by_InputBuffer on empty buffer');

# reset
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
$ib->next();

$c->annotate_InputBuffer($ib);
my $vf = $ib->buffer->[0];

is_deeply($vf->{existing}, [
  {
    'gnomAD_ASJ' => 'T:0',
    'phenotype_or_disease' => 0,
    'gnomAD_SAS' => 'T:0',
    'SAS' => 'T:0',
    'gnomAD' => 'T:0.0003478',
    'failed' => 0,
    'gnomAD_AMR' => 'T:0.0003236',
    'AA' => 'T:0.004998',
    'somatic' => 0,
    'AFR' => 'T:0.003',
    'strand' => 1,
    'gnomAD_FIN' => 'T:0',
    'allele_string' => 'C/T',
    'gnomAD_AFR' => 'T:0.004643',
    'minor_allele_freq' => '0.0010',
    'AMR' => 'T:0.0014',
    'gnomAD_OTH' => 'T:0',
    'gnomAD_EAS' => 'T:0',
    'EUR' => 'T:0',
    'clin_sig' => '',
    'EAS' => 'T:0',
    'end' => 25585733,
    'gnomAD_NFE' => 'T:1.886e-05',
    'variation_name' => 'rs142513484',
    'minor_allele' => 'T',
    'EA' => 'T:0',
    'start' => 25585733,
    'pubmed' => '',
    'matched_alleles' => [
      {
        'a_index' => 0,
        'a_allele' => 'T',
        'b_allele' => 'T',
        'b_index' => 0
      }
    ],
  }
], 'annotate_InputBuffer');

is(scalar (grep {$_->{existing}} @{$ib->buffer}), 132, 'annotate_InputBuffer count annotated');

$vf->{start}++;
delete $vf->{existing};
$ib->{buffer} = [$vf];
$c->annotate_InputBuffer($ib);

is($vf->{existing}, undef, 'annotate_InputBuffer - miss by one');

# construct one to test phenotype_or_disease and clin_sig
$ib = get_ib([qw(21 25891796 . C T . . .)]);
$c->annotate_InputBuffer($ib);

is_deeply(
  $ib->buffer->[0]->{existing},
  [
    {
      'gnomAD_ASJ' => 'T:0',
      'phenotype_or_disease' => '1',
      'gnomAD_SAS' => 'T:0',
      'SAS' => 'T:0',
      'gnomAD' => 'T:9.75e-05',
      'failed' => 0,
      'gnomAD_AMR' => 'T:0.0005957',
      'AA' => '',
      'somatic' => 0,
      'AFR' => 'T:0',
      'strand' => 1,
      'gnomAD_FIN' => 'T:0',
      'allele_string' => 'C/T',
      'gnomAD_AFR' => 'T:0',
      'minor_allele_freq' => '0.0002',
      'AMR' => 'T:0.0014',
      'gnomAD_OTH' => 'T:0.0001823',
      'gnomAD_EAS' => 'T:0',
      'EUR' => 'T:0',
      'clin_sig' => 'uncertain_significance,not_provided,pathogenic',
      'EAS' => 'T:0',
      'end' => 25891796,
      'gnomAD_NFE' => 'T:2.687e-05',
      'variation_name' => 'rs63750066',
      'minor_allele' => 'T',
      'EA' => '',
      'start' => 25891796,
      'pubmed' => '1303275,15365148',
      'matched_alleles' => [
        {
          'a_index' => 0,
          'a_allele' => 'T',
          'b_allele' => 'T',
          'b_index' => 0
        }
      ],
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


# test some nastiness
no warnings 'qw';

$ib = get_ib([qw(21 8987005 . A AGCG . . .)]);
$c->annotate_InputBuffer($ib);
is_deeply(
  $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
  [
    {
      'a_index' => 0,
      'a_allele' => 'GCG',
      'b_allele' => 'GCG',
      'b_index' => 0
    }
  ],
  'nastiness 1'
);

$ib = get_ib([qw(21 8987004 . TA C,TAGCG . . .)]);
$c->annotate_InputBuffer($ib);
is_deeply(
  $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
  [
    {
      'a_index' => 1,
      'a_allele' => 'TAGCG',
      'b_allele' => 'GCG',
      'b_index' => 0
    }
  ],
  'nastiness 2'
);

$ib = get_ib([qw(21 8987004 . TAT TAGCGT . . .)]);
$c->annotate_InputBuffer($ib);
is_deeply(
  $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
  [
    {
      'a_index' => 0,
      'a_allele' => 'AGCGT',
      'b_allele' => 'GCG',
      'b_index' => 0
    }
  ],
  'nastiness 3'
);

$ib = get_ib([qw(21 8987004 . TAT TAGCGT,TAGTGT . . .)]);
$c->annotate_InputBuffer($ib);
is_deeply(
  $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
  [
    {
      'a_index' => 0,
      'a_allele' => 'AGCGT',
      'b_allele' => 'GCG',
      'b_index' => 0
    },
    {
      'a_index' => 1,
      'a_allele' => 'AGTGT',
      'b_allele' => 'GTG',
      'b_index' => 1
    }
  ],
  'nastiness 4'
);


# test old_maf setting
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
$ib->next();

$c->{old_maf} = 1;
$c->clean_cache();
$c->annotate_InputBuffer($ib);
is($ib->buffer->[0]->{existing}->[0]->{AMR}, 0.0014, 'old_maf');
$c->{old_maf} = 0;

# done
done_testing();

sub get_ib {
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      valid_chromosomes => [21],
      file => $test_cfg->create_input_file(shift)
    })
  });
  $ib->next;
  return $ib;
}
