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

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File');


SKIP: {
  no warnings 'once';
  
  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Bio::DB::HTS::Tabix module not available', 22 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

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

  # overlap with insertions
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585746 25585745 -/C)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.vcf.gz' => [
        {
          'name' => 'ins2'
        }
      ]
    },
    'annotate_InputBuffer - overlap - insertion'
  );


  # entry with no ID defaults to reporting coords
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585753 25585753 C/T)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;
  $as->report_coords(0);
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.vcf.gz' => [
        {
          'name' => '21:25585753-25585753'
        }
      ]
    },
    'annotate_InputBuffer - overlap - no ID reports coords'
  );

  
  

  # report_coords
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585733 25585733 C/T + rs142513484)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;

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
  $as->fields(['FOO', 'GOO', 'NOVALUE']);
  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'test.vcf.gz' => [
        { name => 'test1', fields => { 'FOO' => 'BAR', 'GOO' => 'CAR,STAR', 'NOVALUE' => 1, 'FILTER' =>'PASS' } }
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
        { name => 'test1', fields => { 'FOO' => 'BAR', 'GOO' => 'CAR,STAR', 'NOVALUE' => 1, 'FILTER' =>'PASS' } }
      ],
      'foo' => [
        { name => 'test1', allele => 'T', fields => { 'GOO' => 'CAR', 'FILTER' =>'PASS' } }
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
        { name => 'test1', allele => 'A', fields => {'GOO' => 'CAR','FILTER' => 'PASS'} },
        { name => 'test1', allele => 'C', fields => {'GOO' => 'STAR','FILTER' => 'PASS'} }
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
        { name => 'del1', allele => '-', fields => {'GOO' => 'B','FILTER' => 'PASS'} },
        { name => 'del2', allele => '-', fields => {'FILTER' => 'SEGDUP,RF'} }
      ]
    },
    'annotate_InputBuffer - deletion'
  );

  # mixed
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([
        [qw(21 25585741 25585740 -/C +)],
        [qw(21 25585740 25585740 A/C +)]
      ]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;

  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 'ins1', allele => 'C', fields => {'GOO' => 'YAR','FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - mixed - insertion'
  );
  is_deeply(
    $ib->buffer->[1]->{_custom_annotations},
    {
      'foo' => [
        { name => 'ins1', allele => 'C', fields => {'GOO' => 'ZAR', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - mixed - snp'
  );

  # test trimming of input
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([
        [qw(21 25585740 25585740 G/GC + trim_1_l)],
        [qw(21 25585739 25585740 AG/AGC + trim_2_l)],
        [qw(21 25585740 25585741 AG/CG + trim_1_r)],
        [qw(21 25585740 25585742 AGT/CGT + trim_2_r)],
        [qw(21 25585738 25585741 GGAG/GGCG + trim_both)],
      ]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;

  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 'ins1', allele => 'GC', fields => {'GOO' => 'YAR','FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - trim input - insertion 1'
  );
  is_deeply(
    $ib->buffer->[1]->{_custom_annotations},
    {
      'foo' => [
        { name => 'ins1', allele => 'AGC', fields => {'GOO' => 'YAR', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - trim input - insertion 2'
  );
  is_deeply(
    $ib->buffer->[2]->{_custom_annotations},
    {
      'foo' => [
        { name => 'ins1', allele => 'CG', fields => {'GOO' => 'ZAR', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - trim input - snp 1'
  );
  is_deeply(
    $ib->buffer->[3]->{_custom_annotations},
    {
      'foo' => [
        { name => 'ins1', allele => 'CGT', fields => {'GOO' => 'ZAR', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - trim input - snp 2'
  );
  is_deeply(
    $ib->buffer->[4]->{_custom_annotations},
    {
      'foo' => [
        { name => 'ins1', allele => 'GGCG', fields => {'GOO' => 'ZAR', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - trim input - snp 3'
  );

  # test trimming of VCF alleles
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([
        [qw(21 25585771 25585772 TT/- +)],
        [qw(21 25585772 25585773 TT/- +)],
      ]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;

  $as->fields(['FOO']);
  $as->annotate_InputBuffer($ib);
  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 'comp1', allele => '-', fields => {'FOO' => 'CAR', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - trim VCF - deletion 1'
  );
  is_deeply(
    $ib->buffer->[1]->{_custom_annotations},
    {
      'foo' => [
        { name => 'comp1', allele => '-', fields => {'FOO' => 'CAR', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - trim VCF - deletion 2'
  );



  # SV deletion should only have overlap type applied
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VCF->new({
      config => $cfg,
      file => $test_cfg->create_input_file([qw(21 25585735 . G . . . SVTYPE=DEL;SVLEN=2)]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 'del2', fields => { 'FILTER' => 'SEGDUP,RF'} },
      ]
    },
    'annotate_InputBuffer - SV reverts to overlap'
  );

  # test mismatches in num of ALTs to num INFO chunks
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
    config => $cfg,
    parser => Bio::EnsEMBL::VEP::Parser::VEP_input->new({
      config => $cfg,
      file => $test_cfg->create_input_file([
        [qw(21 25585760 25585760 C/T +)],
        [qw(21 25585761 25585761 C/T +)]
      ]),
      valid_chromosomes => [21]
    })
  });
  $ib->next;

  $as->fields(['FOO']);
  $as->annotate_InputBuffer($ib);

  is_deeply(
    $ib->buffer->[0]->{_custom_annotations},
    {
      'foo' => [
        { name => 'refinc', allele => 'T', fields => {'FOO' => 0.1, 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - REF included in INFO chunks'
  );
  
  is_deeply(
    $ib->buffer->[1]->{_custom_annotations},
    {
      'foo' => [
        { name => 'countwrong', allele => 'T', fields => {'FOO' => '0.7,0.2,0.1', 'FILTER' => 'PASS'} },
      ]
    },
    'annotate_InputBuffer - INFO chunk count doesnt match ALT count'
  );
}

done_testing();
