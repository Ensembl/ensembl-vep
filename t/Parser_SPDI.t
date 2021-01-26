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
use_ok('Bio::EnsEMBL::VEP::Parser::SPDI');

# need to get a config object and DB connection for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

throws_ok {
  Bio::EnsEMBL::VEP::Parser::SPDI->new({
    config => Bio::EnsEMBL::VEP::Config->new({offline => 1}),
    file => $test_cfg->create_input_file('21:25000869:G:T')
  });
} qr/Cannot use SPDI format in offline mode/, 'throw without DB';

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 3 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens');

  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    warning_file => 'STDERR',
  });

  my $p = Bio::EnsEMBL::VEP::Parser::SPDI->new({
    config => $cfg,
    file => $test_cfg->create_input_file('21:25000869:G:T'),
    valid_chromosomes => [21],
  });

  is(ref($p), 'Bio::EnsEMBL::VEP::Parser::SPDI', 'class ref');

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
    'allele_string' => 'G/T',
    'ancestral_allele' => undef,
    'map_weight' => 1,
    'chr' => '21',
    '_source_id' => undef,
    'analysis' => undef,
    'end' => 25000870,
    'seq_region_end' => 25000870,
    'minor_allele_frequency' => undef,
    'overlap_consequences' => undef,
    'minor_allele' => undef,
    'start' => 25000870,
    'seq_region_start' => 25000870
  }, 'Bio::EnsEMBL::Variation::VariationFeature' );

  $vf = $p->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'spdi genomic');

  my $spdi_1 = Bio::EnsEMBL::VEP::Parser::SPDI->new({
    config => $cfg,
    file => $test_cfg->create_input_file('LRG_485:6673:G:A'),
    valid_chromosomes => [21,'LRG_485'],
  });

  my $spdi_expected = bless( {
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
    'allele_string' => 'G/A',
    'ancestral_allele' => undef,
    'map_weight' => 1,
    'chr' => 'LRG_485',
    '_source_id' => undef,
    'analysis' => undef,
    'end' => 6674,
    'seq_region_end' => 6674,
    'minor_allele_frequency' => undef,
    'overlap_consequences' => undef,
    'minor_allele' => undef,
    'start' => 6674,
    'seq_region_start' => 6674
  }, 'Bio::EnsEMBL::Variation::VariationFeature' );

  my $spdi_vf = $spdi_1->next();
  delete($spdi_vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($spdi_vf, $spdi_expected, 'spdi genomic LRG');

  1;
};


done_testing();
