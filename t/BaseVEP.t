# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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


## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::BaseVEP');

my $bv = Bio::EnsEMBL::VEP::BaseVEP->new();
ok($bv, 'new is defined');

is(ref($bv), 'Bio::EnsEMBL::VEP::BaseVEP', 'check class');

throws_ok { Bio::EnsEMBL::VEP::BaseVEP->new('not a hashref') } qr/Can\'t use .+ as a HASH ref/, 'new with invalid arg';

throws_ok { Bio::EnsEMBL::VEP::BaseVEP->new({config => {}}) } qr/was expected to be .*Bio::EnsEMBL::VEP::Config/, 'new with invalid config object';



## METHOD TESTS
###############

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new();
ok($cfg, 'get new config object');

$bv = Bio::EnsEMBL::VEP::BaseVEP->new({config => $cfg});
ok($bv, 'new with config object');

is(ref($bv->config), 'Bio::EnsEMBL::VEP::Config', 'config method');

is($bv->param('species'), 'homo_sapiens', 'param get');
is($bv->param('species', 'mouse'), 'mouse', 'param set');
throws_ok { $bv->param() } qr/No param/, 'param without key';
$bv->param('species', 'homo_sapiens');

is($bv->species, 'homo_sapiens', 'species get');
is($bv->species('human'), 'human', 'species set');
$bv->species('homo_sapiens');

is_deeply(
  $bv->stats,
  bless({
    stats => {
      counters => {}
    },
    no_stats => undef,
    _config => $bv->config,
  }, 'Bio::EnsEMBL::VEP::Stats'),
  'stats'
);

# get_adaptor should work offline for some var adaptors using new_fake
$bv->param('offline', 1);
$bv->param('database', 0);
is(ref($bv->get_adaptor('variation', 'VariationFeature')), 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor', 'get_adaptor - offline');

# we can also test throws offline
throws_ok { $bv->get_adaptor() } qr/No adaptor group specified/, 'get_adaptor no group';
throws_ok { $bv->get_adaptor('core') } qr/No adaptor type specified/, 'get_adaptor no type';

# add_shortcuts copies params to $self to speed up lookups
$bv = Bio::EnsEMBL::VEP::BaseVEP->new({
  config => Bio::EnsEMBL::VEP::Config->new({
      %$cfg_hash, test => 'hello'
  })
});
$bv->add_shortcuts('test');
is($bv->{test}, 'hello', 'add_shortcuts');

$bv = Bio::EnsEMBL::VEP::BaseVEP->new({
  config => Bio::EnsEMBL::VEP::Config->new({
      %$cfg_hash, test => 'hello'
  })
});
$bv->add_shortcuts(['test']);
is($bv->{test}, 'hello', 'add_shortcuts arrayref');

throws_ok { $bv->add_shortcuts('test') } qr/add_shortcuts would overwrite value/, 'add_shortcuts overwrite';

is($bv->fasta_db, undef, 'fasta_db - no file');

$bv = Bio::EnsEMBL::VEP::BaseVEP->new({config => Bio::EnsEMBL::VEP::Config->new($cfg_hash)});
is($bv->get_slice('1'), undef, 'get_slice - no fasta_db or database');

use_ok('Bio::EnsEMBL::VEP::Runner');
my $runner = Bio::EnsEMBL::VEP::Runner->new({%$cfg_hash, offline => 1, database => 0, input_file => $test_cfg->{test_vcf}});
$runner->get_all_AnnotationSources();
$bv = Bio::EnsEMBL::VEP::BaseVEP->new({config => $runner->config});

my $fasta_db = $bv->fasta_db;
ok(
  ref($fasta_db) eq 'Bio::DB::HTS::Faidx' || ref($fasta_db) eq 'Bio::DB::Fasta',
  'fasta_db'
);

is_deeply(
  $bv->chr_lengths,
  {21 => 46709983,
   MT => 16569 },
  'chr lengths'
);

is_deeply(
  $bv->get_slice('1'),
  bless( {
    'adaptor' => undef,
    'is_fake' => 1,
    'end' => 1,
    'seq_region_name' => '1',
    'coord_system' => bless( {
      'adaptor' => undef,
      'dbID' => undef,
      'top_level' => 0,
      'version' => undef,
      'name' => 'chromosome',
      'default' => 0,
      'sequence_level' => 0,
      'rank' => 1
    }, 'Bio::EnsEMBL::CoordSystem' ),
    'strand' => 1,
    'seq_region_length' => 1,
    'seq' => undef,
    'start' => 1
  }, 'Bio::EnsEMBL::Slice' ),
  'get_slice - fasta_db'
);

is_deeply($bv->chromosome_synonyms, {}, 'chromosome_syonyms - empty');
my $syns = $bv->chromosome_synonyms($test_cfg->{chr_synonyms});
is(ref($syns), 'HASH', 'chromosome_syonyms - ref');
is_deeply(
  $syns->{21}, {
    'CM000683.2' => 1,
    'NC_000021.9' => 1,
    'chr21' => 1
  },
  'chromosome_syonyms - content check'
);

# fasta_dir
$runner = Bio::EnsEMBL::VEP::Runner->new({
  %$cfg_hash,
  offline => 1,
  database => 0,
  input_file => $test_cfg->{test_vcf},
  fasta_dir => $test_cfg->{fasta_dir}
});
$bv = Bio::EnsEMBL::VEP::BaseVEP->new({config => $runner->config});

ok($bv->fasta_db, 'fasta_dir - fasta_db OK');
is($bv->param('fasta'), $test_cfg->{fasta_dir}.'/Homo_sapiens.GRCh38.toplevel.test.fa', 'fasta_dir - fasta param set');


## status_msg tests require we mess with STDOUT
###############################################

# status_msg prints to STDOUT
no warnings 'once';
open(SAVE, ">&STDOUT") or die "Can't save STDOUT\n"; 

close STDOUT;
my $tmp;
open STDOUT, '>', \$tmp;

$bv->status_msg('Hello');
ok($tmp =~ /Hello/, 'status_msg');

open(STDOUT, ">&SAVE") or die "Can't restore STDOUT\n";


## warning_msg tests require we mess with STDERR
################################################

open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
open STDERR, '>', \$tmp;

# test warning_msg
my $warning_file = $test_cfg->create_input_file('tmp');

$bv = Bio::EnsEMBL::VEP::BaseVEP->new({
  config => Bio::EnsEMBL::VEP::Config->new({
      %$cfg_hash, warning_file => $warning_file, no_progress => 1
  })
});

is(ref($bv->warning_fh), 'FileHandle', 'warning_fh');

ok($bv->warning_msg('test warning message'), 'warning_msg');

# we have to close the fh before we can read the contents of the file
$bv->warning_fh->close();

open IN, $warning_file;
my @contents = map {chomp; $_} <IN>;
close IN;
is($contents[0], 'WARNING: test warning message', 'warning_msg in file');

$bv = Bio::EnsEMBL::VEP::BaseVEP->new({
  config => Bio::EnsEMBL::VEP::Config->new({
      %$cfg_hash, warning_file => 'STDERR', no_progress => 1
  })
});

# test STDERR
$bv->warning_msg('Hello');
ok($tmp =~ /Hello/, 'warning_msg to STDERR');

open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

is($bv->get_database_assembly, undef, 'get_database_assembly - no DB');



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
  skip 'No local database configured', 8 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens');
  
  $bv = Bio::EnsEMBL::VEP::BaseVEP->new({
    config => Bio::EnsEMBL::VEP::Config->new({
      %$cfg_hash,
      %$db_cfg,
      database => 1,
      offline => 0,
      species => 'homo_vepiens',
    })
  });

  ok($bv->registry, 'db - registry');

  is(ref($bv->get_adaptor('core', 'slice')), 'Bio::EnsEMBL::DBSQL::SliceAdaptor', 'get_adaptor slice');

  $bv = Bio::EnsEMBL::VEP::BaseVEP->new({
    config => Bio::EnsEMBL::VEP::Config->new({
      %$cfg_hash,
      registry => $test_cfg->registry_file($multi->{conf}->{core}->{dbname}),
      database => 1,
      offline => 0,
      species => 'homo_vepiens',
    })
  });

  ok($bv->registry, 'db - registry from registry file');

  is(ref($bv->get_adaptor('core', 'slice')), 'Bio::EnsEMBL::DBSQL::SliceAdaptor', 'get_adaptor slice - registry file');

  is(ref($bv->get_slice('21')), 'Bio::EnsEMBL::Slice', 'get_slice - database');

  is($bv->get_database_assembly, 'GRCh38', 'get_database_assembly');

  $multi = Bio::EnsEMBL::Test::MultiTestDB->new('mus_muscuvep');
  
  $bv = Bio::EnsEMBL::VEP::BaseVEP->new({
    config => Bio::EnsEMBL::VEP::Config->new({
      %$cfg_hash,
      %$db_cfg,
      database => 1,
      offline => 0,
      species => 'mus_muscuvep',
    })
  });

  is_deeply(
    $bv->get_adaptor('variation', 'VariationFeature'),
    bless( {
      'species' => 'mus_muscuvep'
    }, 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor' ),
    'get_adaptor - species with no var db gets fake'
  );

  is($bv->get_adaptor('variation', 'Variation'), undef, 'get_adaptor - species with no var db no fake available');

  1;
};


## DONE
#######
done_testing();

