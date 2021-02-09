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
my $base_testing_cfg = $test_cfg->base_testing_cfg;

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::Haplo::Parser::VCF');

use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new($base_testing_cfg);
ok($cfg, 'get new config object');

ok(my $p = Bio::EnsEMBL::VEP::Haplo::Parser::VCF->new({file => $test_cfg->{test_vcf}}), 'new');
is(ref($p), 'Bio::EnsEMBL::VEP::Haplo::Parser::VCF', 'ref');

throws_ok {Bio::EnsEMBL::VEP::Haplo::Parser::VCF->new} qr/No input file/, 'no file';

is_deeply(
  $p->samples,
  {
    'HG00096' => bless( {
      'dbID' => -1,
      'name' => 'HG00096',
      'individual' => bless( {
        'dbID' => -1,
        'name' => 'HG00096',
        'type_individual' => 'outbred'
      }, 'Bio::EnsEMBL::Variation::Individual' ),
      'display' => 'UNDISPLAYABLE'
    }, 'Bio::EnsEMBL::Variation::Sample' )
  },
  'samples'
);

no warnings 'qw';
throws_ok {
  Bio::EnsEMBL::VEP::Haplo::Parser::VCF->new({
    file => $test_cfg->create_input_file([
      [qw(#CHROM POS ID REF ALT QUAL FILTER INFO)],
      [qw(21 1 test1 A G . . .)]
    ])
  })->samples
} qr/no sample/, 'no samples';

no warnings 'qw';
$p = Bio::EnsEMBL::VEP::Haplo::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21], delimiter => "\t"});
is(ref($p->parser), 'Bio::EnsEMBL::IO::Parser::VCF4', 'parser ref');

is_deeply(
  $p->next,
  {
    'alleles' => 'C,T',
    'record' => [
      '21', '25585733', 'rs142513484', 'C', 'T', '.', '.', '.', 'GT', '0|0'
    ],
    'ids' => [
      'rs142513484'
    ],
    'chr' => '21',
    'end' => 25585733,
    'start' => '25585733'
  },
  'next - accounts for header fetch'
);

is($p->next->{ids}->[0], 'rs187353664', 'next again');

# warning_msg prints to STDERR 
no warnings 'once';
open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

close STDERR;
my $tmp;
open STDERR, '>', \$tmp;

is($p->validate_chr({'alleles' => 'C,T', 'chr' => '211', 'end' => 25585733, 'start' => '25585733'}), 0, 'validate chr - chromosome not in valid list'); 
ok($tmp =~ /Chromosome 211 not found in annotation sources/, 'validate chr message - chromosome not in valid list');

is($p->validate_chr({'alleles' => 'C,T', 'chr' => '21', 'end' => 25585733, 'start' => 'x'}), 0, 'validate chr - genomic position not valid');
ok($tmp =~ /WARNING: Start x or end 25585733 coordinate invalid/, 'validate chr message - genomic position invalid');

done_testing();
