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
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript');

my $dir = $test_cfg->{cache_dir};

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
ok($cfg, 'get new config object');

my $c = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
  config => $cfg,
  dir => $dir,
  source_type => 'ensembl'
});
ok($c, 'new is defined');


## METHODS
##########

is($c->serializer_type, 'storable', 'serializer_type');
is($c->file_suffix, 'gz', 'file_suffix');

is($c->get_dump_file_name(1, '1-100'), $dir.'/1/1-100.gz', 'get_dump_file_name');
is($c->get_dump_file_name(1, 1, 100), $dir.'/1/1-100.gz', 'get_dump_file_name with end');

throws_ok { $c->get_dump_file_name() } qr/No chromosome/, 'get_dump_file_name no chromosome';
throws_ok { $c->get_dump_file_name(1) } qr/No region/, 'get_dump_file_name no region';

# deserialization
my $obj = $c->deserialize_from_file(
  $c->get_dump_file_name(
    $test_cfg->{cache_chr},
    $test_cfg->{cache_region}
  )
);
is(ref($obj), 'HASH', 'deserialize_from_file ref 1');
is(ref($obj->{$test_cfg->{cache_chr}}), 'ARRAY', 'deserialize_from_file ref 2');
is(ref($obj->{$test_cfg->{cache_chr}}->[0]), 'Bio::EnsEMBL::Transcript', 'deserialize_from_file ref 3');

# processing deserialized object
my $features = $c->deserialized_obj_to_features($obj);
is(ref($features), 'ARRAY', 'deserialized_object_to_features ref 1');
is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'deserialized_object_to_features ref 2');
is(ref($features->[-1]), 'Bio::EnsEMBL::Transcript', 'deserialized_object_to_features ref 3');
is(scalar @$features, 70, 'deserialized_object_to_features count');

# filter_transcript
is($c->filter_transcript($features->[0]), 1, 'filter_transcript pass');

$c->{gencode_basic} = 1;
is($c->filter_transcript($features->[3]), 0, 'filter_transcript fail gencode_basic');
$c->{gencode_basic} = 0;

$c->{source_type} = 'refseq';
is($c->filter_transcript($features->[0]), 0, 'filter_transcript fail all_refseq');

$features->[0]->{_source_cache} = 'RefSeq';
is($c->filter_transcript($features->[0]), 0, 'filter_transcript fail all_refseq merged');
$c->{source_type} = 'ensembl';

# check filter works on deserialized_obj_to_features
$c->{gencode_basic} = 1;
is(scalar @{$c->deserialized_obj_to_features($obj)}, 47, 'deserialized_object_to_features filtered count');
$c->{gencode_basic} = undef;

is(scalar @{$c->merge_features([@$features, @$features])}, 70, 'merge_features count');

# merge_features does some hacky data restoration
# probably only required for use with old buggy cache files
my @tmp = grep {$_->{_gene_symbol} eq 'MRPL39'} @$features;
delete $tmp[0]->{_gene_hgnc_id};
@tmp = @{$c->merge_features(\@tmp)};
is($tmp[0]->{_gene_hgnc_id}, 'HGNC:14027', 'merge_features restores missing _gene_hgnc_id');

$c->{source_type} = 'refseq';
@tmp = grep {$_->{_gene_symbol} eq 'MRPL39'} @$features;
delete $tmp[0]->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);
@tmp = @{$c->merge_features(\@tmp)};
is($tmp[0]->{_gene_symbol}, 'MRPL39', 'merge_features refseq restores missing _gene_symbol');
is($tmp[0]->{_gene_symbol_source}, 'HGNC', 'merge_features refseq restores missing _gene_symbol_source');
is($tmp[0]->{_gene_hgnc_id}, 'HGNC:14027', 'merge_features refseq restores missing _gene_hgnc_id');
$c->{source_type} = 'ensembl';

$features = $c->get_features_by_regions_from_disk([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);
is(ref($features), 'ARRAY', 'get_features_by_regions_from_disk ref 1');
is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'get_features_by_regions_from_disk ref 2');
is($features->[0]->stable_id, 'ENST00000441009', 'get_features_by_regions_from_disk stable_id');

# now we should be able to retrieve the same from memory
$features = $c->get_features_by_regions_from_memory([[$test_cfg->{cache_chr}, $test_cfg->{cache_s}]]);
is(ref($features), 'ARRAY', 'get_features_by_regions_from_memory ref 1');
is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'get_features_by_regions_from_memory ref 2');
is($features->[0]->stable_id, 'ENST00000441009', 'get_features_by_regions_from_memory stable_id');

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

$features = $c->get_all_features_by_InputBuffer($ib);
is(ref($features), 'ARRAY', 'get_all_features_by_InputBuffer ref 1');
is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'get_all_features_by_InputBuffer ref 2');
is(ref($features->[-1]), 'Bio::EnsEMBL::Transcript', 'get_all_features_by_InputBuffer ref 3');
is($features->[0]->stable_id, 'ENST00000441009', 'get_all_features_by_InputBuffer stable_id');
is(scalar @$features, 70, 'get_all_features_by_InputBuffer count');

# do it again to get them from memory
$features = $c->get_all_features_by_InputBuffer($ib);
is($features->[0]->stable_id, 'ENST00000441009', 'get_all_features_by_InputBuffer again');


$ib->next();
is_deeply($c->get_all_features_by_InputBuffer($ib), [], 'get_all_features_by_InputBuffer on empty buffer');

# reset
$p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
$ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
$ib->next();

$c->annotate_InputBuffer($ib);
my $vf = $ib->buffer->[0];
my $tvs = $vf->get_all_TranscriptVariations;

is(scalar @$tvs, 3, 'annotate_InputBuffer - get_all_TranscriptVariations count');

$vf->_finish_annotation;
is($vf->display_consequence, 'missense_variant', 'annotate_InputBuffer - display_consequence');



## SEREAL
#########

SKIP: {

  eval q{ use Sereal; };
  my $can_use_sereal = $@ ? 0 : 1;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 6 unless $can_use_sereal;

  $c = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
    config => $cfg,
    dir => $test_cfg->{sereal_dir},
    serializer_type => 'sereal',
    source_type => 'ensembl'
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
  is(ref($obj->{$test_cfg->{cache_chr}}), 'ARRAY', 'sereal - deserialize_from_file ref 2');
  is(ref($obj->{$test_cfg->{cache_chr}}->[0]), 'Bio::EnsEMBL::Transcript', 'sereal - deserialize_from_file ref 3');

  is($obj->{$test_cfg->{cache_chr}}->[0]->stable_id, 'ENST00000441009', 'sereal - deserialize_from_file stable_id');

  1;
}

# done
done_testing();
