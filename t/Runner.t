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

my $cfg_hash = $test_cfg->base_testing_cfg;
$cfg_hash->{input_data} = '21 25585733 25585733 C/T + rs142513484';

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Runner');

my $runner = Bio::EnsEMBL::VEP::Runner->new($cfg_hash);
ok($runner, 'new is defined');

is(ref($runner), 'Bio::EnsEMBL::VEP::Runner', 'check class');



## METHOD TESTS
###############

is_deeply(
  $runner->get_all_AnnotationSources(),
  [
    bless( {
      '_config' => $runner->config,
      'cache_region_size' => 1000000,
      'dir' => $test_cfg->{cache_dir},
      'serializer_type' => undef,
      'gencode_basic' => undef,
      'source_type' => 'ensembl',
      'all_refseq' => undef,
      'polyphen' => undef,
      'sift' => undef,
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
    }, 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript' )
  ],
  'get_all_AnnotationSources'
);

# setup_db_connection should return silently in offline mode
ok(!$runner->setup_db_connection(), 'setup_db_connection');

is_deeply($runner->get_Parser, bless({
  '_config' => $runner->config,
  'file' => *Bio::EnsEMBL::VEP::Runner::IN,
  'line_number' => 0,
}, 'Bio::EnsEMBL::VEP::Parser::VEP_input' ), 'get_Parser');

is_deeply($runner->get_InputBuffer, bless({
  '_config' => $runner->config,
  'parser' => bless({
    '_config' => $runner->config,
    'file' => *Bio::EnsEMBL::VEP::Runner::IN,
    'line_number' => 0,
  }, 'Bio::EnsEMBL::VEP::Parser::VEP_input' ),
  'buffer_size' => $runner->param('buffer_size'),
}, 'Bio::EnsEMBL::VEP::InputBuffer' ), 'get_InputBuffer');

my $info = $runner->get_output_header_info;
ok($info->{time}, 'get_output_header_info - time');
ok($info->{api_version} =~ /^\d+$/, 'get_output_header_info - api_version');
ok($info->{vep_version} =~ /^\d+$/, 'get_output_header_info - vep_version');
is(ref($info->{input_headers}), 'ARRAY', 'get_output_header_info - input_headers');

is_deeply(
  $info->{version_data}, 
  {
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
  'get_output_header_info - version_data'
);

is_deeply($runner->get_OutputFactory, bless( {
  '_config' => $runner->config,
  'uniprot' => undef,
  'xref_refseq' => undef,
  'numbers' => undef,
  'polyphen_analysis' => 'humvar',
  'pick_allele' => undef,
  'coding_only' => undef,
  'refseq' => undef,
  'no_escape' => undef,
  'total_length' => undef,
  'per_gene' => undef,
  'sift' => undef,
  'appris' => undef,
  'canonical' => undef,
  'symbol' => undef,
  'pick_order' => $runner->param('pick_order'),
  'pick' => undef,
  'no_intergenic' => undef,
  'gene_phenotype' => undef,
  'flag_pick_allele' => undef,
  'process_ref_homs' => undef,
  'maf_esp' => undef,
  'most_severe' => undef,
  'terms' => 'SO',
  'flag_pick_allele_gene' => undef,
  'flag_pick' => undef,
  'summary' => undef,
  'maf_exac' => undef,
  'polyphen' => undef,
  'gmaf' => undef,
  'biotype' => undef,
  'protein' => undef,
  'domains' => undef,
  'pick_allele_gene' => undef,
  'variant_class' => undef,
  'ccds' => undef,
  'hgvs' => undef,
  'merged' => undef,
  'maf_1kg' => undef,
  'tsl' => undef,
  'pubmed' => undef,
  'header_info' => $info,
}, 'Bio::EnsEMBL::VEP::OutputFactory::VEP_output' ), 'get_OutputFactory');

ok($runner->init, 'init');

is(
  $runner->next_output_line, 
  join("\t", qw(
    rs142513484
    21:25585733
    T
    ENSG00000154719
    ENST00000307301
    Transcript
    3_prime_UTR_variant
    1122
    - - - - -
    IMPACT=MODIFIER;STRAND=-1
  )),
  'next_output_line'
);


## post_setup_checks
####################

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  dir => $test_cfg->{cache_root_dir}.'/sereal',
  hgvs => 1,
  offline => 1,
});
throws_ok {$runner->post_setup_checks} qr/Cannot generate HGVS coordinates in offline mode/, 'post_setup_checks - hgvs + offline with no fasta';


$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  dir => $test_cfg->{cache_root_dir}.'/sereal',
  check_ref => 1,
  offline => 1,
});
throws_ok {$runner->post_setup_checks} qr/Cannot check reference sequences/, 'post_setup_checks - check_ref + offline with no fasta';

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  lrg => 1,
  offline => 1,
});
throws_ok {$runner->post_setup_checks} qr/Cannot map to LRGs in offline mode/, 'post_setup_checks - lrg + offline';

## status_msg tests require we mess with STDOUT
###############################################

# status_msg prints to STDOUT
no warnings 'once';
open(SAVE, ">&STDOUT") or die "Can't save STDOUT\n"; 

close STDOUT;
my $tmp;
open STDOUT, '>', \$tmp;

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  dir => $test_cfg->{cache_root_dir}.'/sereal',
  everything => 1,
  offline => 1,
});
is($runner->param('hgvs'), 1, 'post_setup_checks - everything + offline with no fasta disables hgvs - before');
$runner->post_setup_checks();
is($runner->param('hgvs'), 0, 'post_setup_checks - everything + offline with no fasta disables hgvs - after');
ok($tmp =~ /Disabling --hgvs/, 'post_setup_checks - status_msg for above');
$tmp = '';

foreach my $flag(qw(lrg check_sv check_ref hgvs)) {
  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir}.'/sereal',
    $flag => 1,
    cache => 1,
    offline => 0,
    database => 0,
  });
  $runner->post_setup_checks();
  ok($tmp =~ /Database will be accessed when using --$flag/, 'post_setup_checks - info - '.$flag);
  $tmp = '';
}

open(STDOUT, ">&SAVE") or die "Can't restore STDOUT\n";




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
  
  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
  });

  ok($runner->setup_db_connection(), 'db - setup_db_connection');

  # check it is switching species using aliases
  $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'human_vep',
  });
  $runner->setup_db_connection;

  is($runner->species, 'homo_vepiens', 'db - species alias');

  $info = $runner->get_output_header_info;

  is($info->{db_host}, $db_cfg->{host}, 'get_output_header_info - db_host');
  ok($info->{db_name} =~ /homo_vepiens_core.+/, 'get_output_header_info - db_name');
  ok($info->{db_version}, 'get_output_header_info - db_version');
};

done_testing();
