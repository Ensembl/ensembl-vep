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

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');


SKIP: {
  no warnings 'once';
  
  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 15 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

  ## BASIC TESTS
  ##############

  # use test
  use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File::VCF');

  my $file = $test_cfg->{custom_vcf};

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::VCF->new({file => $file});
  ok($as, 'new is defined');

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Config');

  my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
  ok($cfg, 'get new config object');



  ## TESTS WITH INPUT BUFFER
  ##########################

  use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
  my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({
    config => $cfg,
    file => $test_cfg->create_input_file([qw(21 25585733 rs142513484 C T . . .)]),
    valid_chromosomes => [21]
  });
  ok($p, 'get parser object');

  use_ok('Bio::EnsEMBL::VEP::InputBuffer');
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'check class');

  is(ref($ib->next()), 'ARRAY', 'check buffer next');

  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.vcf.gz' => [
        { name => 'test1' }
      ]
    },
    'annotate_InputBuffer - overlap'
  );

  # report_coords
  delete($ib->buffer->[0]->{_custom_annotations});
  $as->report_coords(1);
  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.vcf.gz' => [
        { name => '21:25585733-25585733' }
      ]
    },
    'annotate_InputBuffer - overlap - report_coords'
  );
  $as->report_coords(0);


  # VCF info fields
  delete($ib->buffer->[0]->{_custom_annotations});
  $as->fields(['FOO', 'GOO']);
  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.vcf.gz' => [
        { name => 'test1', fields => { 'FOO' => 'BAR', 'GOO' => 'CAR,STAR' } }
      ]
    },
    'annotate_InputBuffer - fields'
  );

  # exact type
  $as->type('exact');
  $as->short_name('foo');
  $as->fields(['GOO']);
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.vcf.gz' => [
        { name => 'test1', fields => { 'FOO' => 'BAR', 'GOO' => 'CAR,STAR' } }
      ],
      'foo' => [
        { name => 'test1', allele => 'T', fields => { 'GOO' => 'CAR' } }
      ]
    },
    'annotate_InputBuffer - exact, info keyed on allele'
  );

  # reverse strand input
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585733 25585733 G/A/C -)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;

  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 'test1', allele => 'A', fields => {'GOO' => 'CAR'} },
        { name => 'test1', allele => 'C', fields => {'GOO' => 'STAR'} }
      ]
    },
    'annotate_InputBuffer - exact, multiple, rev strand input'
  );

  # deletion
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585736 25585736 G/- +)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;

  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 'del1', allele => '-', fields => {'GOO' => 'B'} },
        { name => 'del2', allele => '-' }
      ]
    },
    'annotate_InputBuffer - deletion'
  );
}

done_testing();
