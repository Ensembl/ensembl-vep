# Copyright [2016-2026] EMBL-European Bioinformatics Institute
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
  skip 'Bio::DB::BigFile module not available', 25 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_BIGWIG;


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


  ## LARGE-SPAN ZOOM SUMMARY TESTS
  ################################
  # For a large reference span (e.g. a multi-Mb SV) summary_stats min/max are read
  # from the bigWig's precomputed zoom reductions instead of walking every base.
  # This asserts min and max are bit-exact vs the per-base walk, that the
  # per-record output is preserved, and that a combination also needing
  # sum/mean/count falls back to the exact per-base path.
  #
  # A subclass shrinks the exact/zoom window sizes so a small committed fixture
  # (t/testdata/custom/test_large.bw, 3 kb of per-base scores) genuinely
  # exercises the decomposition: with zoom_grid=256 the interior is read from
  # the reduction-160 zoom level, capped by the exact_end=512 end windows.
  {
    no warnings 'once';
    @Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig::SmallWindow::ISA =
      ('Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig');
    *Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig::SmallWindow::zoom_grid = sub { 256 };
  }

  my $large_file = $test_cfg->{custom_bigwig_large};

  my $annotate_large_span = sub {
    my ($fast, $stats) = @_;
    my $class = $fast
      ? 'Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig::SmallWindow'
      : 'Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig';
    my $as_l = $class->new({file => $large_file, format => 'bigwig', config => $cfg});
    $as_l->{summary_stats} = [ @$stats ];
    $as_l->short_name('cons');
    $as_l->type('overlap');

    my $vf_l = { chr => 21, start => 1101, end => 3900 };   # span 2800 bp (> zoom_min_span)
    my $ib_l = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg});
    $ib_l->buffer([$vf_l]);

    if ($fast) {
      $as_l->annotate_InputBuffer($ib_l);                   # may take the zoom path
    }
    else {
      # base-class implementation = the exact per-base walk
      Bio::EnsEMBL::VEP::AnnotationSource::File::annotate_InputBuffer($as_l, $ib_l);
    }
    return $vf_l;
  };

  # summary_stats=min,max -> zoom decomposition, bit-exact vs the per-base walk
  my $exact_vf = $annotate_large_span->(0, [qw(min max)]);
  my $zoom_vf  = $annotate_large_span->(1, [qw(min max)]);

  is(
    $zoom_vf->{_custom_annotations_stats}->{cons}->{max},
    $exact_vf->{_custom_annotations_stats}->{cons}->{max},
    'large span - max bit-exact vs per-base walk'
  );
  is(
    $zoom_vf->{_custom_annotations_stats}->{cons}->{min},
    $exact_vf->{_custom_annotations_stats}->{cons}->{min},
    'large span - min bit-exact vs per-base walk'
  );

  # per-record output parity: same first record, and a truncation marker
  is(
    $zoom_vf->{_custom_annotations}->{cons}->[0]->{name},
    $exact_vf->{_custom_annotations}->{cons}->[0]->{name},
    'large span - first per-record annotation matches per-base walk'
  );
  is(
    $zoom_vf->{_custom_annotations}->{cons}->[-1]->{name}, '...',
    'large span - per-record list truncated with ... marker'
  );

  # a combination that also needs sum/mean/count falls back to the exact per-base
  # path, so its output is identical to the base class (no zoom approximation)
  my $fb_exact = $annotate_large_span->(0, [qw(max mean count)]);
  my $fb_zoom  = $annotate_large_span->(1, [qw(max mean count)]);
  is_deeply(
    $fb_zoom->{_custom_annotations_stats}->{cons},
    $fb_exact->{_custom_annotations_stats}->{cons},
    'large span - stats needing sum/mean/count fall back to exact per-base path'
  );

}


done_testing();
