# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
  skip 'No local database configured', 11 unless $can_use_db;

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

  ## structural variant ID test
  $p = Bio::EnsEMBL::VEP::Parser::ID->new({
    config => $cfg,
    file => $test_cfg->create_input_file('nsv183342'),
    valid_chromosomes => [21],
  });

  is(ref($p), 'Bio::EnsEMBL::VEP::Parser::ID', 'class ref');

  my $expected_sv = bless( {
    'is_somatic' => '0',
    'class_attrib_id' => 10,
    'allele_freq' => undef,
    'outer_start' => undef,
    'seq_region_start' => 25002717,
    '_study_id' => 4279,
    'strand' => 1,
    'seq_region_end' => 25002717,
    'class_SO_term' => 'insertion',
    'allele_string' => undef,
    '_line' => [
      'nsv183342'
    ],
    'outer_end' => undef,
    'chr' => '21',
    'inner_end' => undef,
    '_source_id' => 11,
    'allele_count' => undef,
    'end' => 25002717,
    'length' => undef,
    'breakpoint_order' => undef,
    'inner_start' => undef,
    'start' => 25002717,
    '_structural_variation_id' => 56669931  
  }, 'Bio::EnsEMBL::Variation::StructuralVariationFeature' );

  $vf = $p->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name);

  is_deeply($vf, $expected_sv, 'sv input test');

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

  # sv max size limit
  $tmp = '';
  $p = Bio::EnsEMBL::VEP::Parser::ID->new({
    config => $cfg,
    file => $test_cfg->create_input_file('nsv917902'),
    valid_chromosomes => [21],
  });
  $p->next();
  like($tmp, qr/nsv917902.* variant size .* is bigger than --max_sv_size/, 'SV bigger than max_sv_size warning');

  # restore STDERR
  open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

  1;
};




done_testing();
