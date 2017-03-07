# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
  skip 'Bio::DB::HTS::Tabix module not available', 25 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

  ## BASIC TESTS
  ##############

  # use test
  use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File::GFF');

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Runner');

  my $runner = Bio::EnsEMBL::VEP::Runner->new({%{$test_cfg->base_testing_cfg}, input_file => $test_cfg->{test_vcf}});
  ok($runner, 'get new runner object');

  $runner->init();

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({
    file => $test_cfg->{bam_edit_gff},
    config => $runner->config
  });
  ok($as, 'new is defined');

  ok(!$as->bam, 'no bam provided');

  my $bak = $Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript::CAN_USE_HTS;
  $Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript::CAN_USE_HTS = 0;
  throws_ok {$as->bam($test_cfg->{bam_edit_bam})} qr/Cannot add BAM file/, 'bam - no HTS installed';
  $Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript::CAN_USE_HTS = $bak;

  is(ref($as->bam($test_cfg->{bam_edit_bam})), 'Bio::DB::HTS', 'bam - give filename');

  $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({
    file => $test_cfg->{bam_edit_gff},
    bam => $test_cfg->{bam_edit_bam},
    config => $runner->config
  });

  is(ref($as->bam($test_cfg->{bam_edit_bam})), 'Bio::DB::HTS', 'bam - give filename at new()');  

  # requires fasta db
  ok($as->fasta_db, 'setup fasta_db');

  # requires synonyms
  ok($as->chromosome_synonyms($test_cfg->{chr_synonyms}), 'load synonyms');

  my $l = int(
    $as->fasta_db->length($as->get_source_chr_name(21, 'fasta', [$as->fasta_db->get_all_sequence_ids])) /
    $runner->param('cache_region_size')
  );
  is($l, 46, 'get chrom length');

  my @trs;
  push @trs,
    grep {$_}
    map {$as->lazy_load_transcript($_)}
    @{$as->get_features_by_regions_uncached([['NC_000021.9', $_]])}
    for 1..$l;

  is(scalar @trs, 35, 'get transcripts');

  # test apply edit
  ok($as->apply_edits($trs[0]), 'apply_edits');

  is(
    $trs[0]->spliced_seq,
    'GGGGCGGGGCCGGGCGGCCGTTCAAGGCAGGGGGCGGGGCGTCTCCGAGCGGCGGGGCCAAGGGAGGGCACAACAGCTGCTACCTGAACAGTTTCTGACCCAACAGTTACCCAGCGCCGGACTCGCTGCGCCCCGGCGGCTCTAGGGACCCCCGGCGCCTACACTTAGCTCCGCGCCCGAGAGAATGTTGGACCGACGACACAAGACCTCAGACTTGTGTTATTCTAGCAGCTGAACACACCCCAGGCTCTTCTGACCGGCAGTGGCTCTGGAAGCAGTCTGGTGTATAGAGTTATGGATTCACTACCAGATTCTACTGTATGCTCTTGACAACTATGACCACAATGGTCCACCCACAAATGAATTATCAGGAGTGAACCCAGAGGCACGTATGAATGAAAGTCCTGATCCGACTGACCTGGCGGGAGTCATCATTGAGCTCGGCCCCAATGACAGTCCACAGACAAGTGAATTTAAAGGAGCAACCGAGGAGGCACCTGCGAAAGAAAGCCCACACACAAGTGAATTTAAAGGAGCAGCCCGGGTGTCACCTATCAGTGAAAGTGTGTTAGCACGACTTTCCAAGTTTGAAGTTGAAGATGCTGAAAATGTTGCTTCATATGACAGCAAGATTAAGAAAATTGTGCATTCAATTGTATCATCCTTTGCATTTGGACTATTTGGAGTTTTCCTGGTCTTACTGGATGTCACTCTCATCCTTGCCGACCTAATTTTCACTGACAGCAAACTTTATATTCCTTTGGAGTATCGTTCTATTTCTCTAGCTATTGCCTTATTTTTTCTCATGGATGTTCTTCTTCGAGTATTTGTAGAAAGGAGACAGCAGTATTTTTCTGACTTATTTAACATTTTAGATACTGCCATTATTGTGATTCTTCTGCTGGTTGATGTCGTTTACATTTTTTTTGACATTAAGTTGCTTAGGAATATTCCCAGATGGACACATTTACTTCGACTTCTACGACTTATTATTCTGTTAAGAATTTTTCATCTGTTTCATCAAAAAAGACAACTTGAAAAGCTGATAAGAAGGCGGGTTTCAGAAAACAAAAGGCGATACACAAGGGATGGATTTGACCTAGACCTCACTTACGTTACAGAACGTATTATTGCTATGTCATTTCCATCTTCTGGAAGGCAGTCTTTCTATAGAAATCCAATCAAGGAAGTTGTGCGGTTTCTAGATAAGAAACACCGAAACCACTATCGAGTCTACAATCTATGCAGTGAAAGAGCTTACGATCCTAAGCACTTCCATAATAGGGTCGTTAGAATCATGATTGATGATCATAATGTCCCCACTCTACATCAGATGGTGGTTTTCACCAAGGAAGTAAATGAGTGGATGGCTCAAGATCTTGAAAACATCGTAGCGATTCACTGTAAAGGAGGCACAGATAGAACAGGAACTATGGTTTGTGCCTTCCTTATTGCCTCTGAAATATGTTCAACTGCAAAGGAAAGCCTGTATTATTTTGGAGAAAGGCGAACAGATAAAACCCACAGCGAAAAATTTCAGGGAGTAGAAACTCCTTCTCAGAAGAGATATGTTGCATATTTTGCACAAGTGAAACATCTCTACAACTGGAATCTCCCTCCAAGACGGATACTCTTTATAAAACACTTCATTATTTATTCGATTCCTCGTTATGTACGTGATCTAAAAATCCAAATAGAAATGGAGAAAAAGGTTGTCTTTTCCACTATTTCATTAGGAAAATGTTCGGTACTTGATAACATTACAACAGACAAAATATTAATTGATGTATTCGACGGTCTACCTCTGTATGATGATGTGAAAGTGCAGTTTTTCTATTCGAATCTTCCTACATACTATGACAATTGCTCATTTTACTTCTGGTTGCACACATCTTTTATTGAAAATAACAGGCTTTATCTACCAAAAAATGAATTGGATAATCTACATAAACAAAAAGCACGGAGAATTTATCCATCAGATTTTGCCGTGGAGATACTTTTTGGCGAGAAAATGACTTCCAGTGATGTTGTAGCTGGATCCGATTAAGTATAGCTCCCCCTTCCCCTTCTGGGAAAGAATTATGTTCTTTCCAACCCTGCCACATGTTCATATATCCTAAATCTATCCTAAATGTTCCTTGAAGTATTTATTTATGTTTATATATGTTTATATATGTTCTTCATAAATCTATTACATATATATAGATAAAATGTCAGTGTCTTTCTTCTTTTTTGAAGGTACATACTTCACAGTTTCCATCTTGTATTTACCTAAATTTGGAACACAGTATGCAGGAACATCAGGCATCATTTTGAAGAACTTTGAAATAGAACTTTCAAAGCAAAAACTTGAGATATTCAAAAATTGTGATTCCCCAACACCCACTTACTTACATATTTTTGGTCAAGTATTTATTCCCTGCTTTGGTTCACATTTGTAATTTCAGATTTATTAGAAAGCTAATTTATATATTTTTCTTCCCCTTCATTTACAAACTGTCTGTTAACAGATTGGCAACCAAGACTAAGTTTTAAATCCAGAAGAGAGAAATGTTTCACGCAAGCAGTCCCCCAACTCCCAACACACACTCCCTTTCATTCCGAGCACTAAATAGGTAGCTTACTTGACAAGCTTCCTTTAATTACACACAGACATGGAGGTTGGGGTAAGAAGATGGTGGTAAATATGAAGATAAGTAATCTTTAATAACTTCTGCTTTTGTATAAAATTGTAAGTGAAGTGAAAAGAAATATTCTTAGAGTAA',
    'apply_edits - spliced_seq'
  );

  is_deeply(
    $trs[0]->get_all_Attributes,
    [
      bless( {
        'value' => '256 256 A',
        'name' => 'RNA Edit',
        'description' => 'Edit from bam_edit.bam, op=X len=1 mapped to 256-256',
        'code' => '_rna_edit'
      }, 'Bio::EnsEMBL::Attribute' ),
      bless( {
        'value' => '1547 1547 G',
        'name' => 'RNA Edit',
        'description' => 'Edit from bam_edit.bam, op=X len=1 mapped to 1547-1547',
        'code' => '_rna_edit'
      }, 'Bio::EnsEMBL::Attribute' )
    ],
    'apply_edits - edit attribs'
  );

  is($trs[0]->{_bam_edit_status}, 'ok', 'apply_edits - _bam_edit_status');

  my ($tr) = grep {$_->stable_id eq 'NM_001286476.1'} @trs;
  ok($tr, 'NM_001286476.1 - get');
  ok($as->apply_edits($tr), 'NM_001286476.1 - apply_edits');

  is_deeply(
    $tr->get_all_Attributes,
    [
      bless( {
        'value' => '2966 2965 AAAAAAAAAAAAA',
        'name' => 'RNA Edit',
        'description' => 'Edit from bam_edit.bam, op=I len=13 mapped to 2966-2965',
        'code' => '_rna_edit'
      }, 'Bio::EnsEMBL::Attribute' ),
      bless( {
        'value' => '2941 2940 A',
        'name' => 'RNA Edit',
        'description' => 'Edit from bam_edit.bam, op=I len=1 mapped to 2941-2940',
        'code' => '_rna_edit'
      }, 'Bio::EnsEMBL::Attribute' ),
      bless( {
        'value' => '2618 2617 G',
        'name' => 'RNA Edit',
        'description' => 'Edit from bam_edit.bam, op=I len=1 mapped to 2618-2617',
        'code' => '_rna_edit'
      }, 'Bio::EnsEMBL::Attribute' ),
      bless( {
        'value' => '920 919 C',
        'name' => 'RNA Edit',
        'description' => 'Edit from bam_edit.bam, op=I len=1 mapped to 920-919',
        'code' => '_rna_edit'
      }, 'Bio::EnsEMBL::Attribute' )
    ],
    'NM_001286476.1 - apply_edits - edit attribs'
  );

  $as->apply_edits($_) for @trs;
  is((scalar grep {$_->{_bam_edit_status} eq 'ok'} @trs), 35, 'apply_edits - check all');

  ($tr) =
    grep {$_ && $_->stable_id eq 'NM_001286476.1'}
    map {$as->lazy_load_transcript($_)}
    @{$as->get_features_by_regions_uncached([['NC_000021.9', $_]])}
    for 1..$l;

  ok($tr, 'NM_001286476.1 - reload to force fail');

  substr($tr->get_all_Exons->[0]->{_seq_cache}, 1, 1, 'N');
  ok($as->apply_edits($tr), 'NM_001286476.1 - failed apply_edits returns ok');
  is($tr->{_bam_edit_status}, 'failed', 'NM_001286476.1 - _bam_edit_status failed');
  is(scalar @{$tr->get_all_Attributes}, 0, 'NM_001286476.1 - no edit attribs');
}

done_testing();
