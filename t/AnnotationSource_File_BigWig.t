# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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

use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');

SKIP: {

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  no warnings 'once';
  skip 'Bio::DB::BigFile module not available', 20 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_BIGWIG;


  ## BASIC TESTS
  ##############

  # use test
  use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig');

  my $file = $test_cfg->{custom_bigwig};

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Config');

  my $cfg = Bio::EnsEMBL::VEP::Config->new($test_cfg->base_testing_cfg);
  ok($cfg, 'get new config object');

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig->new({file => $file, format => 'bigwig', config => $cfg});
  ok($as, 'new is defined');
  

  throws_ok {Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig->new({file => 'foo', format => 'bigwig', config => $cfg})->parser} qr/Failed to open/, 'new with invalid file throws';


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
  $as->{summary_stats} = [ 'min', 'mean', 'max', 'sum', 'count' ];

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.bw' => [
        { name => 10 },
      ]
    },
    'annotate_InputBuffer - overlap'
  );

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations_stats},
    {},
    'annotate_InputBuffer - no summary statistics'
  );

  # exact type
  $as->type('exact');
  $as->short_name('foo');
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {

      'test.bw' => [
        { name => 10 },
      ],
      'foo' => [
        { name => 10, score => 10 },
      ]
    },
    'annotate_InputBuffer - exact, additive'
  );

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations_stats},
    {
      'foo' => {
        'min'   => 10,
        'max'   => 10,
        'mean'  => 10,
        'sum'   => 10,
        'count' => 1,
      },
    },
    'annotate_InputBuffer - single score, summary statistics'
  );

  # out by one
  delete($ib->buffer->[0]->{_custom_annotations});
  
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585733 . . <DEL> . . END=25592852)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next();
  
  # overlap multiple scores
  $as->type('overlap');
  $as->short_name('foo');
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 20, score => 20 },
        { name => 30, score => 30 },
        { name => 11, score => 11 },
        { name => 21, score => 21 },
      ]
    },
    'annotate_InputBuffer - overlap, multiple scores'
  );

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations_stats},
    {
      'foo' => {
        'min'   => 11,
        'max'   => 30,
        'mean'  => 20.5,
        'sum'   => 82,
        'count' => 4,
      },
    },
    'annotate_InputBuffer - multiple scores, summary statistics'
  );

  delete($ib->buffer->[0]->{_custom_annotations});

  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585732 rs142513484 C T . . .)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next();

  $as->annotate_InputBuffer($ib);
  ok(!$ib->buffer->[0]->{_custom_annotations}, 'annotate_InputBuffer - out by 1 (5\')');



  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585735 rs142513484 C T . . .)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next();

  $as->annotate_InputBuffer($ib);
  ok(!$ib->buffer->[0]->{_custom_annotations}, 'annotate_InputBuffer - out by 1 (3\')');



  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25592821 rs142513484 C T . . .)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next();

  $as->type('overlap');
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        {
          'name' => '11',
          'score' => '11',
        }
      ]
    },
    'overlap fixedStep'
  );

  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25592821 rs142513484 C T . . .)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next();

  $as->type('overlap');
  $as->report_coords(1);
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        {
          'name' => '21:25592820-25592822',
          'score' => '11',
        }
      ]
    },
    'get scores even when reporting coords'
  );

}


done_testing();
