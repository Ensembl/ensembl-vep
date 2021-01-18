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

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 24 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens');
  
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
    'ancestral_allele' => undef,
    'map_weight' => 1,
    'chr' => '21',
    '_source_id' => undef,
    'analysis' => undef,
    'end' => 25585733,
    'seq_region_end' => 25585733,
    'minor_allele_frequency' => undef,
    'overlap_consequences' => undef,
    'minor_allele' => undef,
    'start' => 25585733,
    'seq_region_start' => 25585733
  }, 'Bio::EnsEMBL::Variation::VariationFeature' );

  $vf = $p->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'genomic');


  # the following return -ve strand equivalent
  $expected->{strand} = -1;
  $expected->{allele_string} = 'G/A';

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENST00000352957.8:c.991G>A'),
    valid_chromosomes => [21],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'coding');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENST00000307301.11:c.*18G>A'),
    valid_chromosomes => [21],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'coding - UTR');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENSP00000284967.6:p.Ala331Thr'),
    valid_chromosomes => [21],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'protein');

  $p = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENST00000352957.8:c.991N>A'),
    valid_chromosomes => [21],
  });
  $p->{lookup_ref} = 1;
  $vf = $p->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'coding - lookup_ref');

  # multiple
  my @vfs;
  $p = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENSP00000284967.6:p.Glu2Asp'),
    valid_chromosomes => [21],
  });
  $p->{ambiguous_hgvs} = 1;

  while($vf = $p->next) {
    push @vfs, $vf;
  }
  is_deeply(
    [sort map {$_->{allele_string}} @vfs],
    ['G/C', 'G/T'],
    'protein - multiple'
  );


  # capture warning
  no warnings 'once';
  open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

  close STDERR;
  open STDERR, '>', \$tmp;

  # multiple from gene input
  @vfs = ();
  $p = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('JAM2:c.721A>T'),
    valid_chromosomes => [21],
  });
  $p->{ambiguous_hgvs} = 1;

  while($vf = $p->next) {
    push @vfs, $vf;
  }
  is_deeply(
    [sort {$a <=> $b} map {$_->{start}} @vfs],
    [25706002, 25706002, 25709954, 25712347],
    'gene - multiple'
  );
  ok($tmp =~ /invalid use of gene or protein identifier/, 'using gene/protein identifier warns');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('21:k.25585733C>T'),
    valid_chromosomes => [21],
  })->next();
  is($vf, undef, 'invalid HGVS type');
  ok($tmp =~ /Unable to parse HGVS notation/, 'invalid HGVS type warning msg');

  $tmp = '';
  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file("21:k.25585733C>T\nENSP00000284967.6:p.Ala331Thr"),
    valid_chromosomes => [21],
  })->next();
  is(ref($vf), 'Bio::EnsEMBL::Variation::VariationFeature', 'skip past invalid HGVS type');
  ok($tmp =~ /Unable to parse HGVS notation/, 'invalid HGVS type warning msg');

  $tmp = '';
  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENSP00000284967.6:p.Glu2Asp'),
    valid_chromosomes => [21],
  })->next();
  is($vf, undef, 'ambiguous protein without ambiguous_hgvs');
  ok($tmp =~ /Could not uniquely determine nucleotide change/, 'ambiguous protein without ambiguous_hgvs msg');

  open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";


  ## REFSEQ
  #########

  # should fetch RefSeq OK by switching core adaptor type automatically
  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('NM_017446.3:c.991G>A'),
    valid_chromosomes => [21],
  })->next();

  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'refseq coding');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('NP_059142.2:p.Ala331Thr'),
    valid_chromosomes => [21],
  })->next();

  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'refseq protein');


  # now try with refseq db as primary
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

  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'refseq primary - refseq coding');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('ENST00000352957.8:c.991G>A'),
    valid_chromosomes => [21],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'refseq primary - ENST coding');




  ## LRG
  ######

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
  });

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

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('LRG_485:g.6674G>A'),
    valid_chromosomes => [21, 'LRG_485'],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'LRG genomic');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('LRG_485t1:c.121G>A'),
    valid_chromosomes => [21, 'LRG_485'],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'LRG coding');

  $vf = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('LRG_485p1:p.Val41Met'),
    valid_chromosomes => [21, 'LRG_485'],
  })->next();
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'LRG protein');


  # test mapping LRG<->chr

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    lrg => 1,
  });

  $p = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('LRG_485:g.6674G>A'),
    valid_chromosomes => [21, 'LRG_485'],
  });

  my $expected_chr = bless( {
    'source' => undef,
    'is_somatic' => undef,
    'display' => undef,
    'clinical_significance' => undef,
    'dbID' => undef,
    'minor_allele_count' => undef,
    '_variation_id' => undef,
    'evidence' => undef,
    'seqname' => undef,
    'strand' => -1,
    'class_SO_term' => undef,
    'map_weight' => 1,
    'allele_string' => 'G/A',
    'ancestral_allele' => undef,
    '_source_id' => undef,
    'chr' => '21',
    'end' => 43774705,
    'seq_region_end' => 43774705,
    'analysis' => undef,
    'minor_allele_frequency' => undef,
    'overlap_consequences' => undef,
    'minor_allele' => undef,
    'start' => 43774705,
    'seq_region_start' => 43774705
  }, 'Bio::EnsEMBL::Variation::VariationFeature' );

  $vf = $p->next;
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected, 'LRG mapping 1');

  $vf = $p->next;
  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply($vf, $expected_chr, 'LRG mapping 2');


  ## check species with no var db, should use fake adaptor(s)
  $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens_coreonly');
  
  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens_coreonly',
    warning_file => 'STDERR',
  });

  $p = Bio::EnsEMBL::VEP::Parser::HGVS->new({
    config => $cfg,
    file => $test_cfg->create_input_file('21:g.25585733C>T '),
    valid_chromosomes => [21],
  });

  $vf = $p->next();

  is_deeply(
    $p->get_adaptor('variation', 'VariationFeature'),
    bless( {
      'species' => 'homo_vepiens_coreonly'
    }, 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor' ),
    'core db only - check adaptor'
  );

  delete($vf->{$_}) for qw(adaptor variation slice variation_name _line);
  is_deeply(
    $vf,
    bless( {
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
      'ancestral_allele' => undef,
      'map_weight' => 1,
      'chr' => '21',
      '_source_id' => undef,
      'analysis' => undef,
      'end' => 25585733,
      'seq_region_end' => 25585733,
      'minor_allele_frequency' => undef,
      'overlap_consequences' => undef,
      'minor_allele' => undef,
      'start' => 25585733,
      'seq_region_start' => 25585733
    }, 'Bio::EnsEMBL::Variation::VariationFeature' ),
    'core db only - genomic hgvs'
  );
  
  1;
};




done_testing();
