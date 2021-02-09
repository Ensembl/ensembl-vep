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
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat');

my $dir = $test_cfg->{cache_dir};

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
ok($cfg, 'get new config object');

my $c = Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat->new({
  config => $cfg,
  dir => $dir,
  cache_region_size => 1000000,
});
ok($c, 'new is defined');


## METHODS
##########

is($c->serializer_type, 'storable', 'serializer_type');
is($c->file_suffix, 'gz', 'file_suffix');

is($c->get_dump_file_name(1, '1-100'), $dir.'/1/1-100_reg.gz', 'get_dump_file_name');
is($c->get_dump_file_name(1, 1, 100), $dir.'/1/1-100_reg.gz', 'get_dump_file_name with end');

throws_ok { $c->get_dump_file_name() } qr/No chromosome/, 'get_dump_file_name no chromosome';
throws_ok { $c->get_dump_file_name(1) } qr/No region/, 'get_dump_file_name no region';

is_deeply($c->get_available_cell_types, [], 'get_available_cell_types - empty');

$c->{available_cell_types} = ['foo', 'bar'];
$c->{cell_type} = ['foo'];
ok($c->check_cell_types, 'check_cell_types');

$c->{cell_type} = ['boo'];
throws_ok {$c->check_cell_types} qr/Cell type .* unavailable/, 'check_cell_types - fail';

# deserialization
my $obj = $c->deserialize_from_file(
  $c->get_dump_file_name(
    $test_cfg->{cache_chr},
    $test_cfg->{cache_region}
  )
);

is(ref($obj), 'HASH', 'deserialize_from_file ref 1');
is(ref($obj->{$test_cfg->{cache_chr}}), 'HASH', 'deserialize_from_file ref 2');

is_deeply(
  [sort keys %{$obj->{$test_cfg->{cache_chr}}}],
  ['MotifFeature', 'RegulatoryFeature'],
  'deserialize_from_file keys'
);

is(
  ref($obj->{$test_cfg->{cache_chr}}->{RegulatoryFeature}->[0]),
  'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
  'deserialize_from_file ref 3'
);
is(
  ref($obj->{$test_cfg->{cache_chr}}->{MotifFeature}->[0]),
  'Bio::EnsEMBL::Funcgen::MotifFeature',
  'deserialize_from_file ref 4'
);

# processing deserialized object
my $features = $c->deserialized_obj_to_features($obj);
is(ref($features), 'ARRAY', 'deserialized_object_to_features ref 1');
is(ref($features->[0]), 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', 'deserialized_object_to_features ref 2');
is(ref($features->[-1]), 'Bio::EnsEMBL::Funcgen::MotifFeature', 'deserialized_object_to_features ref 3');
is(scalar @$features, 1771, 'deserialized_object_to_features count');

is(scalar @{$c->merge_features([@$features, @$features])}, 1605, 'merge_features count');

$features = $c->get_features_by_regions_uncached([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);

is(ref($features), 'ARRAY', 'get_features_by_regions_from_disk ref 1');
is(ref($features->[0]), 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', 'get_features_by_regions_from_disk ref 2');
is($features->[0]->stable_id, 'ENSR00000660046', 'get_features_by_regions_from_disk stable_id');


# now we should be able to retrieve the same from memory
$features = $c->get_features_by_regions_cached([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);
is(ref($features), 'ARRAY', 'get_features_by_regions_from_memory ref 1');
is(ref($features->[0]), 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', 'get_features_by_regions_from_memory ref 2');
is($features->[0]->stable_id, 'ENSR00000660046', 'get_features_by_regions_from_memory stable_id');

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
is(ref($features->[0]), 'Bio::EnsEMBL::Funcgen::RegulatoryFeature', 'get_all_features_by_InputBuffer ref 2');
is(ref($features->[-1]), 'Bio::EnsEMBL::Funcgen::MotifFeature', 'get_all_features_by_InputBuffer ref 3');
is($features->[0]->stable_id, 'ENSR00001054995', 'get_all_features_by_InputBuffer stable_id');
is(scalar @$features, 532, 'get_all_features_by_InputBuffer count');

$features = $c->get_all_features_by_InputBuffer($ib);
is($features->[0]->stable_id, 'ENSR00001054995', 'get_all_features_by_InputBuffer again');

$ib->next();
is_deeply($c->get_all_features_by_InputBuffer($ib), [], 'get_all_features_by_InputBuffer on empty buffer');

# reset
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
$ib->next();

$c->annotate_InputBuffer($ib);
my ($vf) = grep {$_->variation_name eq 'rs3989369'} @{$ib->buffer};
my $rfvs = $vf->get_all_RegulatoryFeatureVariations;

is(scalar @$rfvs, 1, 'annotate_InputBuffer - get_all_RegulatoryFeatureVariations count');

$vf->_finish_annotation;
is($vf->display_consequence, 'regulatory_region_variant', 'annotate_InputBuffer - display_consequence');



## SEREAL
#########

SKIP: {

  eval q{ use Sereal; };
  my $can_use_sereal = $@ ? 0 : 1;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Sereal not installed', 8 unless $can_use_sereal;

  $c = Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat->new({
    config => $cfg,
    dir => $test_cfg->{sereal_dir},
    serializer_type => 'sereal',
  });

  is($c->serializer_type, 'sereal', 'sereal - serializer_type');

  is($c->file_suffix, 'sereal', 'file_suffix');

  # deserialization
  my $obj = $c->deserialize_from_file(
    $c->get_dump_file_name(
      $test_cfg->{cache_chr},
      $test_cfg->{cache_region}
    )
  );

  is(ref($obj), 'HASH', 'sereal - deserialize_from_file ref 1');
  is(ref($obj->{$test_cfg->{cache_chr}}), 'HASH', 'sereal - deserialize_from_file ref 2');

  is_deeply(
    [sort keys %{$obj->{$test_cfg->{cache_chr}}}],
    ['MotifFeature', 'RegulatoryFeature'],
    'sereal - deserialize_from_file keys'
  );

  is(
    ref($obj->{$test_cfg->{cache_chr}}->{RegulatoryFeature}->[0]),
    'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
    'sereal - deserialize_from_file ref 3'
  );

  is(
    ref($obj->{$test_cfg->{cache_chr}}->{MotifFeature}->[0]),
    'Bio::EnsEMBL::Funcgen::MotifFeature',
    'sereal - deserialize_from_file ref 4'
  );

  is(
    $obj->{$test_cfg->{cache_chr}}->{RegulatoryFeature}->[0]->stable_id,
    'ENSR00001565774',
    'sereal - deserialize_from_file stable_id'
  );

  1;
}



# done
done_testing();

