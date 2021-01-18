# Copyright [2016-2021] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
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

my ($vf, $tmp, $expected);

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Parser::ID');

# need to get a config object and DB connection for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

throws_ok {
  Bio::EnsEMBL::VEP::Parser::ID->new({
    config => Bio::EnsEMBL::VEP::Config->new({offline => 1}),
    file => $test_cfg->create_input_file('rs123')
  });
} qr/Cannot use ID format in offline mode/, 'throw without DB';

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 4 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens') if $can_use_db;
  
  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
  });

  my $p = Bio::EnsEMBL::VEP::Parser::ID->new({
    config => $cfg,
    file => $test_cfg->create_input_file('rs142513484'),
    valid_chromosomes => [21],
  });

  is(ref($p), 'Bio::EnsEMBL::VEP::Parser::ID', 'class ref');

  $expected = bless( {
    'is_somatic' => '0',
    'display' => '1',
    'clinical_significance' => [],
    'minor_allele_count' => '5',
    'evidence' => [
      'Multiple_observations',
      'Frequency',
      '1000Genomes',
      'ESP',
      'ExAC'
    ],
    '_variation_id' => '28751744',
    'strand' => '1',
    'class_SO_term' => 'SNV',
    'allele_string' => 'C/T',
    'ancestral_allele' => 'C',
    'map_weight' => '1',
    'chr' => '21',
    '_source_id' => '1',
    'end' => 25585733,
    'seq_region_end' => 25585733,
    'flank_match' => '1',
    'minor_allele_frequency' => '0.000998403',
    'minor_allele' => 'T',
    'start' => 25585733,
    'seq_region_start' => 25585733,
    '_line' => ['rs142513484'],
    '_source_name' => 'dbSNP',
  }, 'Bio::EnsEMBL::Variation::VariationFeature' );

  $vf = $p->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name);

  is_deeply($vf, $expected, 'basic input test');


  my $tmp;
  no warnings 'once';
  open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

  close STDERR;
  open STDERR, '>', \$tmp;

  # missing ID
  $vf = Bio::EnsEMBL::VEP::Parser::ID->new({
    config => $cfg,
    file => $test_cfg->create_input_file('rs699'),
  })->next;
  is($vf, undef, 'missing ID');

  ok($tmp =~ /No variant found with ID/, 'missing ID warning msg');

  # skip past missing ID
  $tmp = '';
  $vf = Bio::EnsEMBL::VEP::Parser::ID->new({
    config => $cfg,
    file => $test_cfg->create_input_file("rs699\nrs142513484"),
  })->next;
  is($vf->{variation_name}, 'rs142513484', 'skip past missing ID');

  ok($tmp =~ /No variant found with ID/, 'missing ID warning msg');

  my $file = $test_cfg->create_input_file();
  open OUT, ">$file";
  print OUT "\n \nrs142513484 \n";
  close OUT;

  $p = Bio::EnsEMBL::VEP::Parser::ID->new({
    config => $cfg,
    file => $file,,
    valid_chromosomes => [21],
  });

  $vf = $p->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name);

  is_deeply($vf, $expected, 'parse input file with empty lines');

  ok($tmp =~ /Skipped 2 empty lines/, 'empty line warning');

  # restore STDERR
  open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

  1;
};




done_testing();
