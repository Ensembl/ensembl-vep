# Copyright [2016] EMBL-European Bioinformatics Institute
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
  'chr' => 21,
  'phenotype_or_disease' => 0,
  'ExAC_AMR' => 'T:0.000173',
  'SAS' => 'T:0.0000',
  'failed' => 0,
  'ExAC_NFE' => 'T:0',
  'AA' => 'T:0.005',
  'somatic' => 0,
  'ExAC_SAS' => 'T:0',
  'AFR' => 'T:0.0030',
  'strand' => 1,
  'allele_string' => 'C/T',
  'ExAC_Adj' => 'T:0.0004133',
  'minor_allele_freq' => '0.0010',
  'ExAC_FIN' => 'T:0',
  'AMR' => 'T:0.0014',
  'EUR' => 'T:0.0000',
  'clin_sig' => undef,
  'EAS' => 'T:0.0000',
  'end' => 25585733,
  'ExAC' => 'T:4.119e-04',
  'ExAC_OTH' => 'T:0',
  'ExAC_AFR' => 'T:0.004681',
  'variation_name' => 'rs142513484',
  'ExAC_EAS' => 'T:0',
  'minor_allele' => 'T',
  'EA' => 'T:0',
  'start' => 25585733,
  'pubmed' => undef,
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
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, valid_chromosomes => [21], file => $test_cfg->create_input_file([qw(21 25891796 . C T . . .)])});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next;

  $c->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{existing},
    [
      {
        'phenotype_or_disease' => '1',
        'ExAC_AMR' => 'T:0.0003457',
        'SAS' => 'T:0.0000',
        'failed' => 0,
        'ExAC_NFE' => 'T:2.998e-05',
        'AA' => undef,
        'somatic' => 0,
        'ExAC_SAS' => 'T:0',
        'AFR' => 'T:0.0000',
        'strand' => 1,
        'allele_string' => 'C/T',
        'ExAC_Adj' => 'T:5.768e-05',
        'minor_allele_freq' => '0.0002',
        'ExAC_FIN' => 'T:0',
        'AMR' => 'T:0.0014',
        'chr' => '21',
        'EUR' => 'T:0.0000',
        'clin_sig' => 'not_provided,pathogenic',
        'EAS' => 'T:0.0000',
        'end' => 25891796,
        'ExAC' => 'T:5.765e-05',
        'ExAC_OTH' => 'T:0.001101',
        'ExAC_AFR' => 'T:0',
        'variation_name' => 'rs63750066',
        'ExAC_EAS' => 'T:0',
        'minor_allele' => 'T',
        'EA' => undef,
        'start' => 25891796,
        'pubmed' => undef
      },
      {
        'phenotype_or_disease' => '1',
        'ExAC_AMR' => undef,
        'SAS' => undef,
        'failed' => 0,
        'ExAC_NFE' => undef,
        'AA' => undef,
        'somatic' => 0,
        'ExAC_SAS' => undef,
        'AFR' => undef,
        'strand' => 1,
        'allele_string' => 'HGMD_MUTATION',
        'ExAC_Adj' => undef,
        'minor_allele_freq' => undef,
        'ExAC_FIN' => undef,
        'AMR' => undef,
        'chr' => '21',
        'EUR' => undef,
        'clin_sig' => undef,
        'EAS' => undef,
        'end' => 25891796,
        'ExAC' => undef,
        'ExAC_OTH' => undef,
        'ExAC_AFR' => undef,
        'variation_name' => 'CM930033',
        'ExAC_EAS' => undef,
        'minor_allele' => undef,
        'EA' => undef,
        'start' => 25891796,
        'pubmed' => undef
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
}



SKIP: {

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 3 unless $Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix::CAN_USE_TABIX_PM;

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();
  $vf = $ib->buffer->[0];

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
  delete($c->{_chr_name_map});
  $c->chromosome_synonyms($test_cfg->{chr_synonyms});
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
}

# done
done_testing();
