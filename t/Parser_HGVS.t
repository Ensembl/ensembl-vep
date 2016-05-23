# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use_ok('Bio::EnsEMBL::VEP::Parser::HGVS');

# need to get a config object and DB connection for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

throws_ok {
  Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => Bio::EnsEMBL::VEP::Config->new({offline => 1}),
    file => $test_cfg->create_input_file('21:g.25585733C>T')
  });
} qr/Cannot use HGVS format in offline mode/, 'throw without DB';

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;
  my $can_use_db = $db_cfg && scalar keys %$db_cfg;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 7 unless $can_use_db;

  my $multi;

  if($can_use_db) {
    eval q{
      use Bio::EnsEMBL::Test::TestUtils;
      use Bio::EnsEMBL::Test::MultiTestDB;
      1;
    };

    $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens');
  }
  
  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    warning_file => 'STDERR',
  });

  my $p = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('21:g.25585733C>T'),
    valid_chromosomes => [21],
  });

  is(ref($p), 'Bio::EnsEMBL::VEP::Parser::HGVS', 'class ref');

  $expected = bless( {
    'source' => undef,
    'is_somatic' => undef,
    'clinical_significance' => undef,
    'display' => undef,
    'dbID' => undef,
    'minor_allele_count' => undef,
    'seqname' => undef,
    'strand' => 1,
    'evidence' => undef,
    '_variation_id' => undef,
    'class_SO_term' => undef,
    'allele_string' => 'C/T',
    'map_weight' => 1,
    'chr' => '21',
    '_source_id' => undef,
    'analysis' => undef,
    'end' => 25585733,
    'minor_allele_frequency' => undef,
    'overlap_consequences' => undef,
    'minor_allele' => undef,
    'start' => 25585733
  }, 'Bio::EnsEMBL::Variation::VariationFeature' );

  $vf = $p->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name);
  is_deeply($vf, $expected, 'genomic');


  # the following return -ve strand equivalent
  $expected->{strand} = -1;
  $expected->{allele_string} = 'G/A';

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENST00000352957.8:c.991G>A'),
    valid_chromosomes => [21],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name);
  is_deeply($vf, $expected, 'coding');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENST00000307301.11:c.*18G>A'),
    valid_chromosomes => [21],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name);
  is_deeply($vf, $expected, 'coding - UTR');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENSP00000284967.6:p.Ala331Thr'),
    valid_chromosomes => [21],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name);
  is_deeply($vf, $expected, 'protein');

  # capture warning
  no warnings 'once';
  open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

  close STDERR;
  open STDERR, '>', \$tmp;

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('21:k.25585733C>T'),
    valid_chromosomes => [21],
  })->next();
  is($vf, undef, 'invalid HGVS type');
  ok($tmp =~ /Unable to parse HGVS notation/, 'invalid HGVS type warning msg');

  open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";

  ## REFSEQ
  #########

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    refseq => 1,
  });

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('NM_017446.3:c.991G>A'),
    valid_chromosomes => [21],
  })->next();

  delete($vf->{$_}) for qw(adaptor variation slice variation_name);
  is_deeply($vf, $expected, 'refseq coding');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('NP_059142.2:p.Ala331Thr'),
    valid_chromosomes => [21],
  })->next();

  delete($vf->{$_}) for qw(adaptor variation slice variation_name);
  is_deeply($vf, $expected, 'refseq protein');
  
  1;
};




done_testing();