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

# get_adaptor should work offline for some var adaptors using new_fake
$bv->param('offline', 1);
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
my $warning_file = $test_cfg->create_input_file();

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



## DATABASE TESTS
#################

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;
  my $can_use_db = $db_cfg && scalar keys %$db_cfg;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 4 unless $can_use_db;

  my $multi;

  if($can_use_db) {
    eval q{
      use Bio::EnsEMBL::Test::TestUtils;
      use Bio::EnsEMBL::Test::MultiTestDB;
      1;
    };

    $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens');
  }
  
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

  1;
};


## DONE
#######
done_testing();

