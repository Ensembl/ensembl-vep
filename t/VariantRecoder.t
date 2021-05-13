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

my $cfg_hash = $test_cfg->base_testing_cfg;

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::VariantRecoder');



SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 22 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens') if $can_use_db;
  
  my $vr = Bio::EnsEMBL::VEP::VariantRecoder->new({%$cfg_hash, %$db_cfg, offline => 0, database => 1, species => 'homo_vepiens'});
  ok($vr, 'new is defined');

  is($vr->param('input_data', 'rs142513484'), 'rs142513484', 'set input_data');

  is(ref($vr), 'Bio::EnsEMBL::VEP::VariantRecoder', 'check class');

  ok($vr->init(), 'init');
  ok($vr->{parser}, 'init sets parser');
  ok($vr->{input_buffer}, 'init sets input_buffer');

  $vr->reset();
  ok(!$vr->{parser}, 'reset deletes parser');
  ok(!$vr->{input_buffer}, 'reset deletes input_buffer');
  ok(!$vr->param('input_data'), 'reset deletes input_data');

  throws_ok {$vr->recode()} qr/No input data/, 'recode - no input';

  foreach my $input(qw(
    ENST00000352957.8:c.991G>A
    NM_017446.3:c.991G>A
    ENSP00000284967.6:p.Ala331Thr
    NP_059142.2:p.Ala331Thr
  )) {
    is_deeply(
      $vr->recode($input),
      [
        {
        "A" =>
        {
          "hgvsp" => [
             "ENSP00000284967.6:p.Ala331Thr",
             "NP_059142.2:p.Ala331Thr",
             "XP_011527953.1:p.Ala289Thr"
          ],
          "hgvsc" => [
             "ENST00000307301.11:c.*18G>A",
             "ENST00000352957.8:c.991G>A",
             "NM_017446.3:c.991G>A",
             "NM_080794.3:c.*18G>A",
             "XM_011529651.1:c.865G>A"
          ],
          "input" => $input,
          "id" => [
             "rs142513484"
          ],
          "hgvsg" => [
             "21:g.25585733C>T"
          ],
          "spdi" => [
             "21:25585732:C:T"
          ]
        }
        }
      ],
      'recode - '.$input
    );
  }

  foreach my $input(qw(
    rs142513484
    21:g.25585733C>T
    21:25585732:C:T
  )) {
    is_deeply(
      $vr->recode($input),
      [
        {
        "T" =>
        {
          "hgvsp" => [
             "ENSP00000284967.6:p.Ala331Thr",
             "NP_059142.2:p.Ala331Thr",
             "XP_011527953.1:p.Ala289Thr"
          ],
          "hgvsc" => [
             "ENST00000307301.11:c.*18G>A",
             "ENST00000352957.8:c.991G>A",
             "NM_017446.3:c.991G>A",
             "NM_080794.3:c.*18G>A",
             "XM_011529651.1:c.865G>A"
          ],
          "input" => $input,
          "id" => [
             "rs142513484"
          ],
          "hgvsg" => [
             "21:g.25585733C>T"
          ],
          "spdi" => [
             "21:25585732:C:T"
          ]
        }
        }
      ],
      'recode - '.$input
    );
  }

  is_deeply(
    [map {@{$_->{T}->{hgvsg}}} @{$vr->recode("JAM2:c.721A>T")}],
    [
       "21:g.25706002A>T",
       "21:g.25709954A>T",
       "21:g.25712347A>T"
    ],
    'recode - from gene gives multiple hgvsg'
  );

  my $bak = $vr->param('fields');
  $vr->param('fields', ['hgvsg']);
  is_deeply(
    $vr->recode("rs142513484"),
    [
      {
      "T" =>
        {
          "input" => "rs142513484",
          "hgvsg" => [
             "21:g.25585733C>T"
          ]
        }
      }
    ],
    'recode - limit fields'
  );
  $vr->param('fields', $bak);

  # Test output VCF format
  my $vr_2 = Bio::EnsEMBL::VEP::VariantRecoder->new({%$cfg_hash, %$db_cfg, offline => 0, database => 1, species => 'homo_vepiens', fields => 'spdi,vcf_string'});
  is_deeply(
    $vr_2->recode("rs142513484"),
    [
      {
      "T" =>
        {
          "input" => "rs142513484",
          "vcf_string" => [
             "21-25585733-C-T"
          ],
          "spdi" => [
             "21:25585732:C:T"
          ]
        }
      }
    ],
    'recode - output vcf_string' 
  );

  my $vr_3 = Bio::EnsEMBL::VEP::VariantRecoder->new({%$cfg_hash, %$db_cfg, offline => 0, database => 1, species => 'homo_vepiens', fields => 'spdi,vcf_string'});
  is_deeply(
    $vr_3->recode("ENST00000352957:c.971_973del"),
    [
      {
      "-" =>
        {
          "input" => "ENST00000352957:c.971_973del",
          "vcf_string" => [
             "21-25585750-GTTA-G"
          ],
          "spdi" => [
             "21:25585750:TTA:"
          ]
        }
      }
    ],
    'recode - input HGVS and output vcf_string'
  );

  my $vr_4 = Bio::EnsEMBL::VEP::VariantRecoder->new({%$cfg_hash, %$db_cfg, offline => 0, database => 1, species => 'homo_vepiens', fields => 'spdi,vcf_string,id,hgvsg'});
  is_deeply(
    $vr_4->recode("rs1444184259"),
    [
      {
      "CCCCCCC" =>
        {
          "input" => "rs1444184259",
          "vcf_string" => [
             "21-25639356-C-CC"
          ],
          "spdi" => [
             "21:25639355:CCCCCC:CCCCCCC"
          ],
          "id" => [
             "rs1444184259"
          ],
          "hgvsg" => [
             "21:g.25639361_25639362insC"
          ]
        }
      }
    ],
    'recode - insertion'
  );

  my $vr_5 = Bio::EnsEMBL::VEP::VariantRecoder->new({%$cfg_hash, %$db_cfg, offline => 0, database => 1, species => 'homo_vepiens', fields => 'spdi,vcf_string,hgvsg'});
  is_deeply(
    $vr_5->recode("21:g.25639361_25639362insC"),
    [
      {
      "C" =>
        {
          "input" => "21:g.25639361_25639362insC",
          "vcf_string" => [
             "21-25639361-C-CC"
          ],
          "spdi" => [
             "21:25639361::C"
          ],
          "hgvsg" => [
             "21:g.25639361dup"
          ]
        }
      }
    ],
    'recode - input HGVS genomic insertion'
  );

  my $vr_6 = Bio::EnsEMBL::VEP::VariantRecoder->new({%$cfg_hash, %$db_cfg, offline => 0, database => 1, species => 'homo_vepiens', var_synonyms => 1});
  is_deeply(
    $vr_6->recode("rs142513484"),
    [
      {
      "T" =>
        {
          "input" => "rs142513484",
          "id" => [
             "rs142513484"
          ],
          "spdi" => [
             "21:25585732:C:T"
          ],
          "var_synonyms" => [
             "LSDB: NM_017446.3:c.991G>A"
          ],
          "hgvsp" => [
             "ENSP00000284967.6:p.Ala331Thr",
             "NP_059142.2:p.Ala331Thr",
             "XP_011527953.1:p.Ala289Thr"
          ],
          "hgvsc" => [
             "ENST00000307301.11:c.*18G>A",
             "ENST00000352957.8:c.991G>A",
             "NM_017446.3:c.991G>A",
             "NM_080794.3:c.*18G>A",
             "XM_011529651.1:c.865G>A"
          ],
          "hgvsg" => [
             "21:g.25585733C>T"
          ]
        }
      }
    ],
    'recode - output variant synonyms'
  );

  };

done_testing();
