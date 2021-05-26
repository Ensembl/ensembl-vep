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
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

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
  use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File::GTF');

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Runner');

  my $runner = Bio::EnsEMBL::VEP::Runner->new({%{$test_cfg->base_testing_cfg}, input_file => $test_cfg->{test_vcf}});
  ok($runner, 'get new runner object');

  $runner->init();

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GTF->new({file => $test_cfg->{custom_gtf}, config => $runner->config});
  ok($as, 'new is defined');


  ok($as->fasta_db, 'setup fasta_db');



  ## METHOD TESTS
  ###############

  # _get_records_by_coords
  my $records = $as->_get_records_by_coords(21, 25585733, 25585733);
  is(scalar @$records, 77, '_get_records_by_coords - count');

  is_deeply(
    $records->[0],
    {
      'source' => 'ensembl_havana',
      'chr' => '21',
      'end' => '25585754',
      'phase' => undef,
      'strand' => '-1',
      'type' => 'exon',
      'md5' => 'e5e3f730a46421a4a0d57887c856c42a',
      'attributes' => {
        'transcript_name' => 'MRPL39-002',
        'gene_source' => 'ensembl_havana',
        'transcript_support_level' => '5',
        'gene_version' => '13',
        'exon_number' => '11',
        'havana_gene' => 'OTTHUMG00000078371',
        'havana_transcript_version' => '2',
        'exon_id' => 'ENSE00003528074',
        'tag' => 'basic',
        'transcript_id' => 'ENST00000307301',
        'gene_id' => 'ENSG00000154719',
        'transcript_biotype' => 'protein_coding',
        'transcript_version' => '11',
        'exon_version' => '1',
        'transcript_source' => 'ensembl_havana',
        'gene_biotype' => 'protein_coding',
        'havana_transcript' => 'OTTHUMT00000171195',
        'havana_gene_version' => '4',
        'ccds_id' => 'CCDS33522',
        'gene_name' => 'MRPL39'
      },
      'start' => '25585656'
    },
    '_get_records_by_coords - first content'
  );

  is(
    scalar (grep {!overlap($_->{start}, $_->{end}, 25585733, 25585733)} @$records),
    59,
    '_get_records_by_coords - non-overlapping sub-features count'
  );

  my $trs = [map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($records)}];

  is(scalar @$trs, 4, "_create_transcripts - count");
  is($trs->[0]->stable_id, "ENST00000307301", "stable_id");
  is($trs->[0]->{_gene_stable_id}, "ENSG00000154719", "_gene_stable_id");
  is($trs->[0]->{_gene_symbol}, "MRPL39", "_gene_symbol");
  is($trs->[0]->{_protein}, "ENSP00000305682", "_protein");
  is(scalar @{$trs->[0]->get_all_Exons}, 11, "count exons");
  my $tr_seq = 'ATGGAGGCGCTGGCCATGGGTTCCCGGGCGCTGCGGCTCTGGCTGGTCGCACCCGGTGGCGGGATCAAATGGAGATTTATAGCAACATCGTCAGCTTCTCAGCTGTCACCGACAGAATTGACAGAAATGCGGAATGATCTCTTTAATAAAGAGAAAGCCAGGCAGTTATCATTAACTCCCCGAACTGAGAAGATAGAAGTTAAGCATGTTGGGAAAACTGACCCCGGTACTGTCTTCGTGATGAATAAAAACATTTCAACTCCCTACAGTTGTGCCATGCATTTAAGCGAGTGGTATTGCAGGAAGTCCATTCTGGCTCTGGTGGATGGACAGCCTTGGGACATGTATAAGCCTTTAACAAAGTCCTGTGAAATTAAATTTCTTACTTTCAAAGATTGTGATCCAGGAGAAGTGAATAAGGCATATTGGCGTTCCTGTGCTATGATGATGGGCTGTGTGATAGAGAGGGCATTCAAAGATGAATATATGGTCAATTTGGTCAGAGCTCCAGAAGTTCCTGTAATTTCTGGTGCCTTCTGTTATGACGTAGTTTTGGATAGCAAACTTGATGAGTGGATGCCAACAAAAGAGAACTTACGTTCCTTCACAAAAGATGCTCATGCTTTAATTTATAAAGATCTTCCATTTGAAACTCTGGAAGTTGAAGCAAAAGTGGCATTGGAAATATTTCAACACAGCAAGTACAAAGTAGATTTCATAGAAGAGAAGGCATCTCAGAACCCTGAGAGAATAGTCAAGCTACACAGAATAGGTGACTTCATTGATGTGAGTGAGGGCCCTCTTATTCCAAGAACAAGTATTTGTTTCCAGTATGAAGTATCAGCAGTTCACAATCTTCAACCCACCCAGCCAAGTCTCATACGAAGATTCCAGGGCGTGTCTTTACCTGTTCACTTAAGAGCACATTTTACAATATGGGATAAGCTATTGGAAAGATCTCGGAAAATGACTCCATTTCCCATTCTCCTTCTATTTACTACGCAGTCATTCTTCACTACCTCGCCTGAGTCGTACCTCCTCCATGGAACAGTCTCAGAGTAA';
  is($trs->[0]->translateable_seq, $tr_seq, "translateable_seq");
  my $prot_seq = 'MEALAMGSRALRLWLVAPGGGIKWRFIATSSASQLSPTELTEMRNDLFNKEKARQLSLTPRTEKIEVKHVGKTDPGTVFVMNKNISTPYSCAMHLSEWYCRKSILALVDGQPWDMYKPLTKSCEIKFLTFKDCDPGEVNKAYWRSCAMMMGCVIERAFKDEYMVNLVRAPEVPVISGAFCYDVVLDSKLDEWMPTKENLRSFTKDAHALIYKDLPFETLEVEAKVALEIFQHSKYKVDFIEEKASQNPERIVKLHRIGDFIDVSEGPLIPRTSICFQYEVSAVHNLQPTQPSLIRRFQGVSLPVHLRAHFTIWDKLLERSRKMTPFPILLLFTTQSFFTTSPESYLLHGTVSE';
  is($trs->[0]->translation->seq, $prot_seq, "translate");

  $trs = [map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($as->_get_records_by_coords(21, 255e5, 26e6))}];
  is_deeply(
    {map {$_->stable_id => $_->biotype} @$trs},
    {
      'ENST00000440126' => 'protein_coding',
      'ENST00000457143' => 'protein_coding',
      'ENST00000415997' => 'protein_coding',
      'ENST00000284971' => 'protein_coding',
      'ENST00000384075' => 'snRNA',
      'ENST00000460679' => 'nonsense_mediated_decay',
      'ENST00000400094' => 'protein_coding',
      'ENST00000600590' => 'antisense',
      'ENST00000464867' => 'retained_intron',
      'ENST00000455275' => 'antisense',
      'ENST00000400099' => 'protein_coding',
      'ENST00000400087' => 'protein_coding',
      'ENST00000307301' => 'protein_coding',
      'ENST00000599572' => 'antisense',
      'ENST00000617755' => 'antisense',
      'ENST00000352957' => 'protein_coding',
      'ENST00000597894' => 'antisense',
      'ENST00000487266' => 'processed_transcript',
      'ENST00000462267' => 'retained_intron',
      'ENST00000357903' => 'protein_coding',
      'ENST00000466453' => 'processed_transcript',
      'ENST00000486002' => 'retained_intron',
      'ENST00000354192' => 'protein_coding',
      'ENST00000456904' => 'antisense',
      'ENST00000419694' => 'lincRNA',
      'ENST00000419219' => 'protein_coding',
      'ENST00000419219_1' => 'protein_coding', # Fake transcript to test exons having 2 parents
      'ENST00000400532' => 'protein_coding',
      'ENST00000596669' => 'retained_intron',
      'ENST00000477351' => 'processed_transcript',
      'ENST00000609365' => 'lincRNA',
      'ENST00000359726' => 'protein_coding',
      'ENST00000358918' => 'protein_coding',
      'ENST00000567517' => 'antisense',
      'ENST00000456917' => 'lincRNA',
      'ENST00000548570' => 'processed_transcript',
      'ENST00000492962' => 'retained_intron',
      'ENST00000471689' => 'retained_intron',
      'ENST00000439274' => 'protein_coding',
      'ENST00000448850' => 'protein_coding',
      'ENST00000420965' => 'processed_pseudogene',
      'ENST00000463070' => 'processed_transcript',
      'ENST00000400093' => 'protein_coding',
      'ENST00000491395' => 'processed_transcript',
      'ENST00000474136' => 'processed_transcript',
      'ENST00000400090' => 'protein_coding',
      'ENST00000400075' => 'protein_coding',
      'ENST00000450769' => 'processed_pseudogene',
      'ENST00000516163' => 'snRNA',
      'ENST00000596385' => 'antisense',
      'ENST00000354828' => 'protein_coding',
      'ENST00000312957' => 'protein_coding',
      'ENST00000346798' => 'protein_coding',
      'ENST00000436405' => 'processed_pseudogene',
      'ENST00000608591' => 'lincRNA',
      'ENST00000480456' => 'protein_coding',
      'ENST00000385060' => 'miRNA',
      'ENST00000348990' => 'protein_coding'
    },
    '_create_transcripts - big fetch check biotypes'
  );


  ## TESTS WITH INPUT BUFFER
  ##########################

  use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
  my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({
    config => $runner->config,
    file => $test_cfg->create_input_file([qw(21 25585733 rs142513484 C T . . .)]),
    valid_chromosomes => [21]
  });
  ok($p, 'get parser object');

  use_ok('Bio::EnsEMBL::VEP::InputBuffer');
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $runner->config, parser => $p});
  is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'check class');

  is(ref($ib->next()), 'ARRAY', 'check buffer next');

  $as->annotate_InputBuffer($ib);
  $ib->finish_annotation();

  is($ib->buffer->[0]->display_consequence, 'missense_variant', 'annotate_InputBuffer - display_consequence');

  ## Test chromosome MT
  my $runner_mt = Bio::EnsEMBL::VEP::Runner->new({%{$test_cfg->base_testing_cfg}, input_file => $test_cfg->{test_vcf_MT}});
  $runner_mt->init();

  my $gtf_mt = Bio::EnsEMBL::VEP::AnnotationSource::File::GTF->new({file => $test_cfg->{custom_gtf_mt}, config => $runner_mt->config});

  my $p_mt = Bio::EnsEMBL::VEP::Parser::VCF->new({
    config => $runner_mt->config,
    file => $test_cfg->create_input_file([qw(MT 4472 var1 T A . . .)]),
    valid_chromosomes => ['MT']
  });

  my $ib_mt = Bio::EnsEMBL::VEP::InputBuffer->new({config => $runner_mt->config, parser => $p_mt});
  is(ref($ib_mt->next()), 'ARRAY', 'check buffer next (MT)');

  $gtf_mt->annotate_InputBuffer($ib_mt);
  $ib_mt->finish_annotation();
  is ($ib_mt->buffer->[0]->get_all_TranscriptVariations->[0]->_codon_table, 2, 'codon table for MT chromosome is correct');
  is ($ib_mt->buffer->[0]->get_all_TranscriptVariations->[0]->pep_allele_string, 'I/M', 'codon table for MT chromosome - check allele');

}


done_testing();
