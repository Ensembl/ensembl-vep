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
      'valid_chromosomes' => [21, 'LRG_485'],
    }, 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript' )
  ],
  'get_all_AnnotationSources'
);

# setup_db_connection should return silently in offline mode
ok(!$runner->setup_db_connection(), 'setup_db_connection');

is_deeply($runner->get_valid_chromosomes, [21, 'LRG_485'], 'get_valid_chromosomes');

is_deeply($runner->get_Parser, bless({
  '_config' => $runner->config,
  'file' => *Bio::EnsEMBL::VEP::Runner::IN,
  'line_number' => 0,
  'check_ref' => undef,
  'chr' => undef,
  'dont_skip' => undef,
  'valid_chromosomes' => {21 => 1, LRG_485 => 1},
  'minimal' => undef,
  'lrg' => undef,
}, 'Bio::EnsEMBL::VEP::Parser::VEP_input' ), 'get_Parser');

is_deeply($runner->get_InputBuffer, bless({
  '_config' => $runner->config,
  'parser' => bless({
    '_config' => $runner->config,
    'file' => *Bio::EnsEMBL::VEP::Runner::IN,
    'line_number' => 0,
    'check_ref' => undef,
    'chr' => undef,
    'dont_skip' => undef,
    'valid_chromosomes' => {21 => 1, LRG_485 => 1},
    'minimal' => undef,
    'lrg' => undef,
  }, 'Bio::EnsEMBL::VEP::Parser::VEP_input' ),
  'buffer_size' => $runner->param('buffer_size'),
  'minimal' => undef,
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
  'plugins' => [],
  'no_stats' => undef,
}, 'Bio::EnsEMBL::VEP::OutputFactory::VEP_output' ), 'get_OutputFactory');

ok($runner->init, 'init');

is_deeply(
  $runner->_buffer_to_output($runner->get_InputBuffer),
  [
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
    join("\t", qw(
      rs142513484
      21:25585733
      T
      ENSG00000154719
      ENST00000352957
      Transcript
      missense_variant
      1033
      991
      331
      A/T
      Gca/Aca
      -
      IMPACT=MODERATE;STRAND=-1
    )),
    join("\t", qw(
      rs142513484
      21:25585733
      T
      ENSG00000260583
      ENST00000567517
      Transcript
      upstream_gene_variant
      -
      -
      -
      -
      -
      -
      IMPACT=MODIFIER;DISTANCE=2407;STRAND=-1
    )),
  ],
  '_buffer_to_output'
);

$runner = Bio::EnsEMBL::VEP::Runner->new($cfg_hash);

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


# plugins
$runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, plugin => ['TestPlugin'], quiet => 1});

is_deeply(
  $runner->get_all_Plugins,
  [
    bless( {
      'params' => [],
      'variant_feature_types' => [
        'VariationFeature'
      ],
      'feature_types' => [
        'Transcript'
      ],
      'version' => '2.3',
      'config' => $runner->config,
    }, 'TestPlugin' )
  ],
  'get_all_Plugins'
);

# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

my $tmp;
close STDERR;
open STDERR, '>', \$tmp;

# plugin doesn't compile
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNoCompile'],
  quiet => 1,
});
ok($runner->get_all_Plugins, 'get_all_Plugins - failed to compile');
ok($tmp =~ /Failed to compile plugin/, 'get_all_Plugins - failed to compile message');

# $runner = Bio::EnsEMBL::VEP::Runner->new({
#   %$cfg_hash,
#   warning_file => 'STDERR',
#   plugin => ['TestPluginNoCompile'],
#   quiet => 1,
#   safe => 1
# });
# throws_ok {$runner->get_all_Plugins} qr/Failed to compile plugin/, 'get_all_Plugins - failed to compile safe die';


# plugin new method fails
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNewFails'],
  quiet => 1,
});
ok($runner->get_all_Plugins, 'get_all_Plugins - new fails');
ok($tmp =~ /Failed to instantiate plugin/, 'get_all_Plugins - new fails message');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNewFails'],
  quiet => 1,
  safe => 1
});
throws_ok {$runner->get_all_Plugins} qr/Failed to instantiate plugin/, 'get_all_Plugins - new fails safe die';


# plugin missing required methods
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNoMethods'],
  quiet => 1,
});
ok($runner->get_all_Plugins, 'get_all_Plugins - missing methods');
ok($tmp =~ /required method/, 'get_all_Plugins - missing methods message');

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  warning_file => 'STDERR',
  plugin => ['TestPluginNoMethods'],
  quiet => 1,
  safe => 1
});
throws_ok {$runner->get_all_Plugins} qr/required method/, 'get_all_Plugins - missing methods safe die';


# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";



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



## FORKING
##########

my $input = [
  [qw(21 25585733 rs142513484 C T . . . GT 0|0)],
  [qw(21 25587701 rs187353664 T C . . . GT 0|0)],
  [qw(21 25587758 rs116645811 G A . . . GT 0|0)],
  [qw(21 25588859 rs199510789 C T . . . GT 0|0)],
];

my $exp = [
  'rs142513484 T ENST00000307301',
  'rs142513484 T ENST00000352957',
  'rs142513484 T ENST00000567517',
  'rs187353664 C ENST00000307301',
  'rs187353664 C ENST00000352957',
  'rs187353664 C ENST00000567517',
  'rs116645811 A ENST00000307301',
  'rs116645811 A ENST00000352957',
  'rs116645811 A ENST00000567517',
  'rs199510789 T ENST00000307301',
  'rs199510789 T ENST00000352957',
  'rs199510789 T ENST00000419219'
];

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->create_input_file($input),
  input_data => undef,
});

my @lines;
while(my $line = $runner->next_output_line) {
  my @split = split("\t", $line);
  push @lines, join(" ", $split[0], $split[2], $split[4]);
}

is_deeply(
  \@lines,
  $exp,
  'fork - check no fork'
);

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->create_input_file($input),
  input_data => undef,
  fork => 2,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  my @split = split("\t", $line);
  push @lines, join(" ", $split[0], $split[2], $split[4]);
}

is_deeply(
  \@lines,
  $exp,
  'fork - check fork (2)'
);

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->create_input_file($input),
  input_data => undef,
  fork => 4,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  my @split = split("\t", $line);
  push @lines, join(" ", $split[0], $split[2], $split[4]);
}

is_deeply(
  \@lines,
  $exp,
  'fork - check fork (4)'
);


# check whole input file
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->{test_vcf},
  input_data => undef,
});

$exp = [];
while(my $line = $runner->next_output_line) {
  push @$exp, $line;
}

# lets check stats aswell
my $exp_stats = $runner->stats->{stats}->{counters};

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->{test_vcf},
  input_data => undef,
  fork => 2,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  push @lines, $line;
}

is_deeply(
  \@lines,
  $exp,
  'fork - full file check (2)'
);

is_deeply(
  $runner->stats->{stats}->{counters},
  $exp_stats,
  'fork - stats (2)'
);

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  input_file => $test_cfg->{test_vcf},
  input_data => undef,
  fork => 4,
});

@lines = ();
while(my $line = $runner->next_output_line) {
  push @lines, $line;
}

is_deeply(
  \@lines,
  $exp,
  'fork - full file check (4)'
);

is_deeply(
  $runner->stats->{stats}->{counters},
  $exp_stats,
  'fork - stats (4)'
);






# warning_msg prints to STDERR
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  fork => 2,
});

$runner->param('warning_file', 'STDERR');

$runner->{_test_warning} = 1;

$runner->next_output_line();
ok($tmp =~ 'TEST WARNING', 'fork - test warning');


$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  fork => 1,
});
$runner->param('warning_file', 'STDERR');
$runner->{_test_die} = 1;

throws_ok {$runner->next_output_line} qr/TEST DIE/, 'fork - test die';

# restore STDERR
open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";



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
