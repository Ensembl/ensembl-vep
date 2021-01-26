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

## BASIC TESTS
##############

use_ok('Bio::EnsEMBL::VEP::AnnotationType::Transcript');
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');


SKIP: {
  no warnings 'once';

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available or Set::IntervalTree not installed', 11 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM && $Bio::EnsEMBL::VEP::AnnotationType::Transcript::CAN_USE_INTERVAL_TREE;

  # use test
  use_ok('Bio::EnsEMBL::VEP::Haplo::AnnotationSource::File::GTF');


  use_ok('Bio::EnsEMBL::VEP::Haplo::Runner');
  my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new({
    %{$test_cfg->base_testing_cfg},
    input_file => $test_cfg->{test_vcf},
    fasta => $test_cfg->{fasta},
    gtf => $test_cfg->{custom_gtf},
    offline => 0,
  });

  ok(
    my $as = Bio::EnsEMBL::VEP::Haplo::AnnotationSource::File::GTF->new({
      file => $test_cfg->{custom_gtf},
      config => $runner->config
    }),
    'new is defined'
  );
  is(ref($as), 'Bio::EnsEMBL::VEP::Haplo::AnnotationSource::File::GTF', 'new ref');

  my $t;
  ok($t = Bio::EnsEMBL::VEP::TranscriptTree->new({annotation_source => $as, config => $runner->config}), 'get TranscriptTree');
  delete($t->{trees});

  # new will have run populate_tree, so should be able to test fetch now
  is_deeply($t->fetch(21, 25585733, 25585733), [], 'fetch before populated empty');

  $as->populate_tree($t);

  is_deeply($t->fetch(21, 25585733, 25585733), [[25585656, 25607517]], 'fetch after populated');

  my $ib = $runner->get_InputBuffer;
  $ib->next();

  is(scalar @{$ib->buffer}, 34, 'input buffer next count');

  my $result = $as->annotate_InputBuffer($ib);

  is(ref($result->[0]), 'Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer', 'annotate_InputBuffer - ref');

  is_deeply(
    [sort map {$_->transcript->stable_id} @$result],
    [qw(ENST00000307301 ENST00000352957  ENST00000419219  ENST00000419219_1)],
    'annotate_InputBuffer - stable_ids'
  );

  is_deeply(
    [
      sort
      map {$_->_hex}
      map {@{$_->get_all_TranscriptHaplotypes}}
      @$result
    ],
    [
      '03f82458a34fda1fdf74803a02cf0a45',
      '03f82458a34fda1fdf74803a02cf0a45',
      '1a040be3cfa8edb6cc091f942d42dcf6',
      '1a040be3cfa8edb6cc091f942d42dcf6',
      '252ea7bfa7d0b130fc696ff25f4fc33c',
      '2ff7fea7ee58c8a776fc97f01f9d5296',
      '41d44949cb17e204c3e76ce6d67aa27d',
      '41d44949cb17e204c3e76ce6d67aa27d',
      '68f1b015c35163ade06017b7b6bce1ab',
      '79d615f301917d96d4cefaea29702508',
      '7aa776c0543c7d064e0d146b60541843',
      '91212e92e2a1940184a2eb77fd115d5c',
      'b3dd19425bfa07a93783baf201f79a70',
      'b3dd19425bfa07a93783baf201f79a70',
      'c6c26cfaa79a2046a684220c04b5b2e4',
      'e5745909c573dfc830946d4e6994a47d'
    ],
    'annotate_InputBuffer - TH hexes'
  );
}

done_testing();
