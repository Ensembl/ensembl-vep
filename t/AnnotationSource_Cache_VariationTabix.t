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

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix');

my $dir = $test_cfg->{sereal_dir};

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
ok($cfg, 'get new config object');

my $c = Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix->new({
  config => $cfg,
  dir => $dir,
  cols => ['chr', @{$test_cfg->{var_cols}}],
});
ok($c, 'new is defined');


## METHODS
##########

is($c->get_dump_file_name(1, '1-100'), $dir.'/1/all_vars.gz', 'get_dump_file_name');

throws_ok { $c->get_dump_file_name() } qr/No chromosome/, 'get_dump_file_name no chromosome';

is($c->delimiter, "\t", 'delimiter');


## TESTS WITH AN INPUT BUFFER
#############################

use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
ok($p, 'get parser object');

my $exp = [{
  'gnomAD_ASJ' => 'T:0',
  'phenotype_or_disease' => 0,
  'gnomAD_SAS' => 'T:0',
  'SAS' => 'T:0',
  'gnomAD' => 'T:0.0003478',
  'failed' => 0,
  'gnomAD_AMR' => 'T:0.0003236',
  'AA' => 'T:0.004998',
  'somatic' => 0,
  'AFR' => 'T:0.003',
  'strand' => 1,
  'gnomAD_FIN' => 'T:0',
  'allele_string' => 'C/T',
  'gnomAD_AFR' => 'T:0.004643',
  'minor_allele_freq' => '0.0010',
  'AMR' => 'T:0.0014',
  'chr' => '21',
  'gnomAD_OTH' => 'T:0',
  'gnomAD_EAS' => 'T:0',
  'EUR' => 'T:0',
  'clin_sig' => undef,
  'EAS' => 'T:0',
  'end' => 25585733,
  'gnomAD_NFE' => 'T:1.886e-05',
  'variation_name' => 'rs142513484',
  'minor_allele' => 'T',
  'EA' => 'T:0',
  'start' => 25585733,
  'pubmed' => undef,
  'matched_alleles' => [
    {
      'a_index' => 0,
      'a_allele' => 'T',
      'b_allele' => 'T',
      'b_index' => 0
    }
  ],
}];

use_ok('Bio::EnsEMBL::VEP::InputBuffer');
my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'check class');

is(ref($ib->next()), 'ARRAY', 'check buffer next');

my ($vf, $vf_hash);

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'tabix binary not available', 3
    unless $Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::CAN_USE_TABIX_CL;

  # the two methods in this class use a hashref of lists of VFs keyed on chr
  $vf = $ib->buffer->[0];
  $vf_hash = {
    $vf->{chr} => [$vf],
  };

  $c->_annotate_cl($vf_hash);

  is_deeply($vf->{existing}, $exp, '_annotate_cl');

  # check synonyms
  $c->chromosome_synonyms($test_cfg->{chr_synonyms});
  $c->{valid_chromosomes} = [21];
  $vf->{chr} = 'NC_000021.9';

  delete $vf->{existing};
  $vf_hash = {
    $vf->{chr} => [$vf],
  };

  $c->_annotate_cl($vf_hash);

  is_deeply($vf->{existing}, $exp, 'chr synonym');

  # no match
  $vf->{chr} = 21;
  $vf->{start}++;

  delete $vf->{existing};
  $vf_hash = {
    $vf->{chr} => [$vf],
  };

  is_deeply($vf->{existing}, undef, 'miss by coord - _annotate_cl');
}

SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module or tabix binary not available', 19
    unless $Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::CAN_USE_TABIX_CL
    or $Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::CAN_USE_TABIX_PM;

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  $c->annotate_InputBuffer($ib);
  $vf = $ib->buffer->[0];

  is_deeply($vf->{existing}, $exp, 'annotate_InputBuffer');

  is(scalar (grep {$_->{existing}} @{$ib->buffer}), 132, 'annotate_InputBuffer count annotated');

  # construct one to test phenotype_or_disease and clin_sig
  $ib = get_ib([qw(21 25891796 . C T . . .)]);
  $c->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{existing},
    [
      {
        'gnomAD_ASJ' => undef,
        'phenotype_or_disease' => '1',
        'gnomAD_SAS' => undef,
        'SAS' => undef,
        'gnomAD' => undef,
        'failed' => 0,
        'gnomAD_AMR' => undef,
        'AA' => undef,
        'somatic' => 0,
        'AFR' => undef,
        'strand' => 1,
        'gnomAD_FIN' => undef,
        'allele_string' => 'HGMD_MUTATION',
        'gnomAD_AFR' => undef,
        'minor_allele_freq' => undef,
        'AMR' => undef,
        'chr' => '21',
        'gnomAD_OTH' => undef,
        'gnomAD_EAS' => undef,
        'EUR' => undef,
        'clin_sig' => undef,
        'EAS' => undef,
        'end' => 25891796,
        'gnomAD_NFE' => undef,
        'variation_name' => 'CM930033',
        'minor_allele' => undef,
        'EA' => undef,
        'start' => 25891796,
        'pubmed' => undef
      },
      {
        'gnomAD_ASJ' => 'T:0',
        'phenotype_or_disease' => '1',
        'gnomAD_SAS' => 'T:0',
        'SAS' => 'T:0',
        'gnomAD' => 'T:9.75e-05',
        'failed' => 0,
        'gnomAD_AMR' => 'T:0.0005957',
        'AA' => undef,
        'somatic' => 0,
        'AFR' => 'T:0',
        'strand' => 1,
        'gnomAD_FIN' => 'T:0',
        'allele_string' => 'C/T',
        'gnomAD_AFR' => 'T:0',
        'minor_allele_freq' => '0.0002',
        'AMR' => 'T:0.0014',
        'chr' => '21',
        'gnomAD_OTH' => 'T:0.0001823',
        'gnomAD_EAS' => 'T:0',
        'EUR' => 'T:0',
        'clin_sig' => 'not_provided,pathogenic',
        'EAS' => 'T:0',
        'end' => 25891796,
        'gnomAD_NFE' => 'T:2.687e-05',
        'variation_name' => 'rs63750066',
        'minor_allele' => 'T',
        'EA' => undef,
        'start' => 25891796,
        'pubmed' => undef,,
        'matched_alleles' => [
          {
            'a_index' => 0,
            'a_allele' => 'T',
            'b_allele' => 'T',
            'b_index' => 0
          }
        ],
      }
    ],
    'annotate_InputBuffer - phenotype_or_disease'
  );


  ## FREQUENCY STUFF
  ##################

  # for now this requires we switch off allele checking
  $c->{no_check_alleles} = 1;

  # new checks
  ok(
    $c->check_frequency_filter,
    'check_frequency_filter'
  );

  my $bak = $c->{freq_pop};
  $c->{freq_pop} = 'foo';
  throws_ok { $c->check_frequency_filter } qr/Invalid population/, 'check_frequency_filter - fails on missing pop';

  $c->{freq_pop} = 'gnomAD_AFR_AF';
  throws_ok { $c->check_frequency_filter } qr/Invalid population/, 'check_frequency_filter - fails on pop with AF added';

  $c->{freq_pop} = 'ExAC_AFR';
  throws_ok { $c->check_frequency_filter } qr/Invalid population/, 'check_frequency_filter - fails on pop for incorrect grouping';

  $c->{freq_pop} = $bak;

  # get_frequency_data
  $vf = $ib->buffer->[0];
  $c->get_frequency_data($vf);

  is_deeply(
    $vf->{_freq_check_freqs},
    {
      '1KG_ALL' => {
        'T' => '0.0002'
      },
    },
    'get_frequency_data - _freq_check_freqs'
  );

  is_deeply(
    $vf->{_freq_check_pass},
    {
      'T' => 0
    },
    'get_frequency_data - _freq_check_pass'
  );

  is($vf->{_freq_check_all_failed}, 1, 'get_frequency_data - _freq_check_all_failed');

  my %orig = %$c;

  # test switching gt_lt
  $c->{freq_gt_lt} = 'lt';
  $c->get_frequency_data($vf);

  is_deeply(
    $vf->{_freq_check_pass},
    {
      'T' => 1
    },
    'get_frequency_data - gt_lt'
  );
  is($vf->{_freq_check_all_passed}, 1, 'get_frequency_data - gt_lt - _freq_check_all_passed');

  $c->{freq_gt_lt} = $orig{freq_gt_lt};

  # test adjusting freq
  $c->{freq_freq} = 0.0001;
  $c->get_frequency_data($vf);

  is_deeply(
    $vf->{_freq_check_pass},
    {
      'T' => 1
    },
    'get_frequency_data - freq'
  );
  is($vf->{_freq_check_all_passed}, 1, 'get_frequency_data - freq - _freq_check_all_passed');
  $c->{freq_freq} = $orig{freq_freq};

  # test another pop
  $c->{freq_pop} = '1KG_AMR';
  delete($vf->{_freq_check_freqs});
  $c->get_frequency_data($vf);

  is_deeply(
    $vf->{_freq_check_freqs},
    {
      '1KG_AMR' => {
        'T' => '0.0014'
      },
    },
    'get_frequency_data - pop _freq_check_freqs'
  );
  $c->{freq_pop} = $orig{freq_pop};

  # test strand switching
  $vf->{strand} = -1;
  $vf->{allele_string} = 'G/A';
  delete($vf->{_freq_check_freqs});
  delete($vf->{_alt_alleles});
  $c->get_frequency_data($vf);

  is_deeply(
    $vf->{_freq_check_freqs},
    {
      '1KG_ALL' => {
        'A' => '0.0002'
      },
    },
    'get_frequency_data - rev strand'
  );


  ## frequency_check_buffer

  # exclude (default)
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, valid_chromosomes => [21], file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next;
  $c->annotate_InputBuffer($ib);
  $c->frequency_check_buffer($ib);
  is(scalar @{$ib->buffer}, 114, 'frequency_check_buffer - exclude count');

  # include
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, valid_chromosomes => [21], file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next;
  $c->annotate_InputBuffer($ib);
  $c->{freq_filter} = 'include';
  $c->frequency_check_buffer($ib);
  is(scalar @{$ib->buffer}, 18, 'frequency_check_buffer - include count');
  $c->{freq_filter} = $orig{freq_filter};

  # test removing single allele
  $c->{freq_freq} = 0.0001;

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, valid_chromosomes => [21], file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next;

  $vf = $ib->buffer->[0];
  $vf->{allele_string} = 'C/T/G';
  $ib->reset_buffer();
  $ib->buffer([$vf]);

  $c->annotate_InputBuffer($ib);
  $c->frequency_check_buffer($ib);
  is($vf->{allele_string}, 'C/G', 'frequency_check_buffer - remove allele 1');

  # include
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, valid_chromosomes => [21], file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next;

  $vf = $ib->buffer->[0];
  $vf->{allele_string} = 'C/T/G';
  $ib->buffer([$vf]);

  $c->{freq_filter} = 'include';
  $c->annotate_InputBuffer($ib);
  $c->frequency_check_buffer($ib);
  is($vf->{allele_string}, 'C/T', 'frequency_check_buffer - remove allele 2');
  $c->{freq_filter} = $orig{freq_filter};

  # reset
  $c->{freq_freq} = $orig{freq_freq};
  $c->{no_check_alleles} = 0;




  # test some nastiness
  no warnings 'qw';

  $ib = get_ib([qw(21 8987005 . A AGCG . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 0,
        'a_allele' => 'GCG',
        'b_allele' => 'GCG',
        'b_index' => 0
      }
    ],
    'nastiness 1'
  );

  $ib = get_ib([qw(21 8987004 . TA C,TAGCG . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 1,
        'a_allele' => 'TAGCG',
        'b_allele' => 'GCG',
        'b_index' => 0
      }
    ],
    'nastiness 2'
  );

  $ib = get_ib([qw(21 8987004 . TAT TAGCGT . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 0,
        'a_allele' => 'AGCGT',
        'b_allele' => 'GCG',
        'b_index' => 0
      }
    ],
    'nastiness 3'
  );

  $ib = get_ib([qw(21 8987004 . TAT TAGCGT,TAGTGT . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 0,
        'a_allele' => 'AGCGT',
        'b_allele' => 'GCG',
        'b_index' => 0
      },
      {
        'a_index' => 1,
        'a_allele' => 'AGTGT',
        'b_allele' => 'GTG',
        'b_index' => 1
      }
    ],
    'nastiness 4'
  );
}



SKIP: {

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 6 unless $Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::CAN_USE_TABIX_PM;

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();
  $vf = $ib->buffer->[0];

  my $tabix_obj = $c->_get_tabix_obj(21);
  is(ref($tabix_obj), 'Bio::DB::HTS::Tabix', '_get_tabix_obj - ref');
  ok($tabix_obj eq $c->_get_tabix_obj(21), '_get_tabix_obj - cache OK');

  delete $c->{_tabix_obj_cache};
  ok($tabix_obj ne $c->_get_tabix_obj(21), '_get_tabix_obj - clear cache new obj');

  $vf_hash = {
    21 => [$vf],
  };
  $c->_annotate_pm($vf_hash);

  is_deeply($vf->{existing}, $exp, '_annotate_pm');


  # check synonyms
  $vf->{chr} = 'NC_000021.9';
  delete $vf->{existing};
  $vf_hash = {
    $vf->{chr} => [$vf],
  };
  $c->chromosome_synonyms($test_cfg->{chr_synonyms});
  $c->{valid_chromosomes} = [21];
  $c->_annotate_pm($vf_hash);

  is_deeply($vf->{existing}, $exp, '_annotate_pm - chr synonym');
  $vf->{chr} = 21;


  $vf->{start}++;
  delete $vf->{existing};
  $vf_hash = {
    21 => [$vf],
  };
  $c->_annotate_pm($vf_hash);

  is_deeply($vf->{existing}, undef, 'miss by coord - _annotate_pm');
}


SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'tabix binary not available', 2
    unless $Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::CAN_USE_TABIX_CL;


  # do the same test with the "opposite" use of perl module vs command line
  $Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::CAN_USE_TABIX_PM = 0;

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  $c->annotate_InputBuffer($ib);
  $vf = $ib->buffer->[0];

  is_deeply($vf->{existing}, $exp, 'annotate_InputBuffer - _annotate_cl');

  is(scalar (grep {$_->{existing}} @{$ib->buffer}), 132, 'annotate_InputBuffer count annotated');


  # test some nastiness
  no warnings 'qw';

  $ib = get_ib([qw(21 8987005 . A AGCG . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 0,
        'a_allele' => 'GCG',
        'b_allele' => 'GCG',
        'b_index' => 0
      }
    ],
    'nastiness 1 - _annotate_cl'
  );

  $ib = get_ib([qw(21 8987004 . TA C,TAGCG . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 1,
        'a_allele' => 'TAGCG',
        'b_allele' => 'GCG',
        'b_index' => 0
      }
    ],
    'nastiness 2 - _annotate_cl'
  );

  $ib = get_ib([qw(21 8987004 . TAT TAGCGT . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 0,
        'a_allele' => 'AGCGT',
        'b_allele' => 'GCG',
        'b_index' => 0
      }
    ],
    'nastiness 3 - _annotate_cl'
  );

  $ib = get_ib([qw(21 8987004 . TAT TAGCGT,TAGTGT . . .)]);
  $c->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{existing}->[0]->{matched_alleles},
    [
      {
        'a_index' => 0,
        'a_allele' => 'AGCGT',
        'b_allele' => 'GCG',
        'b_index' => 0
      },
      {
        'a_index' => 1,
        'a_allele' => 'AGTGT',
        'b_allele' => 'GTG',
        'b_index' => 1
      }
    ],
    'nastiness 4 - _annotate_cl'
  );
}

# done
done_testing();

sub get_ib {
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      valid_chromosomes => [21],
      file => $test_cfg->create_input_file(shift)
    })
  });
  $ib->next;
  return $ib;
}
