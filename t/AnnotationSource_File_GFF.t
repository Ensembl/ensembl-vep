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
  skip 'Bio::DB::HTS::Tabix module not available', 52 unless $Bio::EnsEMBL::VEP::AnnotationSource::File::CAN_USE_TABIX_PM;

  ## BASIC TESTS
  ##############

  # use test
  use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File::GFF');

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Runner');

  my $runner = Bio::EnsEMBL::VEP::Runner->new({%{$test_cfg->base_testing_cfg}, input_file => $test_cfg->{test_vcf}});
  ok($runner, 'get new runner object');

  $runner->init();

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({file => $test_cfg->{custom_gff}, config => $runner->config});
  ok($as, 'new is defined');


  ok($as->fasta_db, 'setup fasta_db');

  throws_ok {
    $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({file => $test_cfg->{custom_gff}, config => $runner->config, type => 'exact'})
  } qr/GXF annotation sources cannot be used as \"exact\" custom annotation type/, 'cant use exact type';



  ## METHOD TESTS
  ###############

  # _get_records_by_coords
  my $records = $as->_get_records_by_coords(21, 25585733, 25585733);
  is(scalar @$records, 75, '_get_records_by_coords - count');

  is_deeply(
    $records->[0],
    {
      'source' => 'ensembl_havana',
      'chr' => '21',
      'end' => '25585754',
      'phase' => undef,
      'strand' => '-1',
      'type' => 'exon',
      'md5' => '4ef31f6bc93e125ab6ed9ed662af743e',
      'attributes' => {
        'ensembl_end_phase' => '-1',
        'version' => '1',
        'exon_id' => 'ENSE00003528074',
        'constitutive' => '0',
        'rank' => '11',
        'Name' => 'ENSE00003528074',
        'ensembl_phase' => '2',
        'Parent' => 'transcript:ENST00000307301'
      },
      'start' => '25585656'
    },
    '_get_records_by_coords - first content'
  );

  is(
    scalar (grep {!overlap($_->{start}, $_->{end}, 25585733, 25585733)} @$records),
    57,
    '_get_records_by_coords - non-overlapping sub-features count'
  );


  my $trs = [map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($records)}];

  is(scalar @$trs, 4, "_create_transcripts - count");
  is($trs->[0]->stable_id, "ENST00000307301", "stable_id");
  is($trs->[0]->{_gene_stable_id}, "ENSG00000154719", "_gene_stable_id");
  is($trs->[0]->{_gene_symbol}, "MRPL39", "_gene_symbol");
  is($trs->[0]->{_protein}, "ENSP00000305682", "_protein");
  is_deeply(
    $trs->[0]->get_all_Attributes,
    [
      bless( {
        'value' => 'tsl5',
        'code' => 'TSL'
      }, 'Bio::EnsEMBL::Attribute' )
    ],
    'attributes'
  );
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


  $records = $as->_get_records_by_coords(21, 25585733, 25585733);
  delete $records->[3]->{attributes}->{gene_id};
  delete $records->[2]->{attributes}->{gene_id};
  is([map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($records)}]->[0]->{_gene_stable_id}, 'ENSG00000154719', 'gene_stable_id from gene object ID not gene_id');

  ## TEST SOME FAILSAFES
  ######################

  # parent/child structure broken
  no warnings 'once';
  open(SAVE, ">&STDERR") or die "Can't save STDERR\n"; 

  my $tmp;
  close STDERR;
  open STDERR, '>', \$tmp;

  $as->_get_parent_child_structure([
    {
      md5 => 1,
      attributes => { parent => 'foo', id => 'test1' }
    }
  ]);
  ok($tmp =~ /Parent entries with the following IDs were not found/, '_get_parent_child_structure - parent not found');

  # invalid child type
  $records = $as->_get_records_by_coords(21, 25585733, 25585733);
  $records->[0]->{type} = 'foo';
  throws_ok { map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($records)} } qr/unexpected type of child record/, '_create_transcripts - unexpected type of child';

  $records = $as->_get_records_by_coords(21, 25585733, 25585733);
  delete $records->[3]->{attributes}->{biotype};
  is(scalar (grep {defined($_)} map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($records)}), 3, 'no biotype skips transcript');
  ok($tmp =~ /Unable to determine biotype/, 'no biotype warning message');

  # overlapping exons
  $records = $as->_get_records_by_coords(21, 25585733, 25585733);
  $records->[0]->{end} += 1e5;
  delete($records->[0]->{attributes}->{rank});
  is(scalar (grep {defined($_)} map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($records)}), 3, 'overlapping exons skips transcript');
  ok($tmp =~ /Failed to add exon to transcript/, 'overlapping exons warning message');

  # missing exons for protein_coding transcript
  my %feature_record = (
     '_children' => [
         { '_parent_id' => [ 'parent_gene_id.1' ],
           'attributes' => { 'Parent' => 'parent_gene_id.1' },
           'chr' => '21',
           'end' => 36705932,
           'phase' => 0,
           'source' => 'test',
           'start' => 36705821,
           'strand' => '-1',
           'type' => 'CDS',
        },
        {  '_parent_id' => [ 'parent_gene_id.1' ],
           'attributes' =>  { 'Parent' => 'parent_gene_id.1'},
           'chr' => '21',
           'end' => 36705714,
           'phase' => 2,
           'source' => 'test',
           'start' => 36705179,
           'strand' => '-1',
           'type' => 'CDS',
        } ],
     '_gene_record' => {
        '_parent_id' => [],
        'attributes' => {
           'ID' => 'parent_gene_id',
         },
        'chr' => '21',
        'end' => 36705932,
        'phase' => undef,
        'source' => 'test',
        'start' => 36705179,
        'strand' => '-1',
        'type' => 'gene',
      },
     '_id' => 'parent_gene_id.1',
     '_parent_id' => [ 'parent_gene_id' ],
     'attributes' => {
        'ID' => 'parent_gene_id.1',
        'Parent' => 'parent_gene_id',
      },
     'chr' => '21',
     'end' => 36705932,
     'phase' => undef,
     'source' => 'test',
     'start' => 36705179,
     'strand' => '-1',
     'type' => 'mRNA',
  );
  my $trans = $as->lazy_load_transcript(\%feature_record, $feature_record{_gene_record});
  ok($tmp =~ /No exons found for protein_coding transcript/, 'no exons warning message');

  # restore STDERR
  open(STDERR, ">&SAVE") or die "Can't restore STDERR\n";



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

  ## transcript filter
  $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({
    file => $test_cfg->{custom_gff},
    config => $runner->config,
    filter => 'stable_id ne ENST00000352957',
  });

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $as->config, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $as->config, parser => $p});
  $ib->next();

  $as->annotate_InputBuffer($ib);
  my $vf = $ib->buffer->[0];
  $vf->_finish_annotation;
  is(scalar (grep {$_->transcript->stable_id eq 'ENST00000352957'} @{$vf->get_all_TranscriptVariations}), 0, 'with filter - filtered transcript absent');
  is($vf->display_consequence, '3_prime_UTR_variant', 'with filter - display_consequence');


  ## get_source_chr_name tests
  ############################

  $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({file => $test_cfg->{custom_gff}, config => $runner->config});

  is($as->get_source_chr_name(21), 21, 'get_source_chr_name - exists, same');
  is($as->get_source_chr_name('chr21'), 21, 'get_source_chr_name - strip chr');

  push @{$as->{valid_chromosomes}}, 'chrFoo';
  is($as->get_source_chr_name('Foo'), 'chrFoo', 'get_source_chr_name - add chr');


  ## REFSEQ FILE
  ##############

  $as = Bio::EnsEMBL::VEP::AnnotationSource::File::GFF->new({file => $test_cfg->{custom_refseq_gff}, config => $runner->config});

  # requires synonyms
  ok($as->chromosome_synonyms($test_cfg->{chr_synonyms}), 'load synonyms');

  # _get_records_by_coords
  $records = $as->_get_records_by_coords('NC_000021.9', 25585733, 25585733);
  is(scalar @$records, 104, '_get_records_by_coords - refseq - count');

  is_deeply(
    $records->[0],
    {
      'source' => 'BestRefSeq',
      'chr' => 'NC_000021.9',
      'end' => '25585754',
      'phase' => undef,
      'strand' => '-1',
      'type' => 'exon',
      'md5' => '71e05169668470578431dfee84361fa7',
      'attributes' => {
        'ID' => 'id1679974',
        'gene' => 'MRPL39',
        'Note' => 'The RefSeq transcript has 1 substitution compared to this genomic sequence',
        'Parent' => 'rna144911',
        'gbkey' => 'mRNA',
        'exception' => 'annotated by transcript or proteomic data',
        'product' => 'mitochondrial ribosomal protein L39%2C transcript variant 1',
        'Dbxref' => 'GeneID:54148,Genbank:NM_017446.3,HGNC:HGNC:14027,HPRD:11371,MIM:611845',
        'transcript_id' => 'NM_017446.3'
      },
      'start' => '25585656'
    },
    '_get_records_by_coords - check first'
  );

  $trs = [map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($records)}];

  is(scalar @$trs, 4, "_create_transcripts - count");
  is($trs->[3]->stable_id, "NM_080794.3", "stable_id");
  is($trs->[3]->{_gene_stable_id}, "54148", "_gene_stable_id");
  is($trs->[3]->{_gene_symbol}, "MRPL39", "_gene_symbol");
  is($trs->[3]->{_protein}, "NP_542984.2", "_protein");
  is(scalar @{$trs->[3]->get_all_Exons}, 11, "count exons");
  is($trs->[3]->translation->seq, $prot_seq, "translate");

  $trs = [map {$as->lazy_load_transcript($_)} @{$as->_create_transcripts($as->_get_records_by_coords('NC_000021.9', 255e5, 26e6))}];

  is_deeply(
    {map {$_->stable_id => $_->biotype} @$trs},
    {
      'XR_937621.2' => 'lncRNA',
      'NR_024092.1' => 'lncRNA',
      'XR_001754985.1' => 'lncRNA',
      'XR_937617.2' => 'lncRNA',
      'NM_021219.3' => 'protein_coding',
      'NM_000484.3' => 'protein_coding',
      'NM_201414.2' => 'protein_coding',
      'XR_937622.2' => 'lncRNA',
      'NM_001003703.1' => 'protein_coding',
      'NM_001320266.1' => 'protein_coding',
      'XR_001755128.1' => 'lncRNA',
      'XR_937537.2' => 'misc_RNA',
      'XM_005260938.4' => 'protein_coding',
      'XM_017028313.1' => 'protein_coding',
      'NR_001458.3' => 'lncRNA',
      'XR_937619.2' => 'lncRNA',
      'XR_001754983.1' => 'lncRNA',
      'NM_001003696.1' => 'protein_coding',
      'XR_001754988.1' => 'lncRNA',
      'NM_080794.3' => 'protein_coding',
      'NM_001270408.1' => 'protein_coding',
      'NM_001136129.2' => 'protein_coding',
      'rna144907' => 'miRNA',
      'NM_001003697.1' => 'protein_coding',
      'NM_201413.2' => 'protein_coding',
      'NM_001136130.2' => 'protein_coding',
      'NM_001136131.2' => 'protein_coding',
      'XR_001754986.1' => 'lncRNA',
      'NM_001204303.1' => 'protein_coding',
      'NM_001204302.1' => 'protein_coding',
      'XR_001755127.1' => 'lncRNA',
      'NM_001270407.1' => 'protein_coding',
      'XM_006724026.3' => 'protein_coding',
      'NM_001685.4' => 'protein_coding',
      'XM_011529520.2' => 'protein_coding',
      'NM_001204301.1' => 'protein_coding',
      'XM_011529651.2' => 'protein_coding',
      'NM_001320267.1' => 'protein_coding',
      'NM_002040.3' => 'protein_coding',
      'NR_030784.1' => 'miRNA',
      'NM_001136016.3' => 'protein_coding',
      'NM_001197297.1' => 'protein_coding',
      'NM_001003701.1' => 'protein_coding',
      'NM_017446.3' => 'protein_coding',
      'XR_001754984.1' => 'lncRNA',
      'XR_001754987.1' => 'lncRNA',
      'rna144908' => 'miRNA',
      'NR_072999.1' => 'misc_RNA',
      'XM_011529521.2' => 'protein_coding'
    },
    '_create_transcripts - big fetch check biotype'
  );

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({
    config => $runner->config,
    file => $test_cfg->create_input_file([qw(21 25585733 rs142513484 C T . . .)]),
    valid_chromosomes => [21]
  });
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $runner->config, parser => $p});
  is(ref($ib->next()), 'ARRAY', 'check buffer next');

  $as->annotate_InputBuffer($ib);
  $ib->finish_annotation();

  is($ib->buffer->[0]->display_consequence, 'missense_variant', 'annotate_InputBuffer - display_consequence');
  # dont_skip test Check that variants that are on seq_regions which are not part of the GFF file are still
  # included in the output if --dont_skip is used
  my $in = qq{21\t25585733\trs142513484\tC\tT\t.\t.\t.
  SEQ_21\t25587758\trs116645811\tG\tA\t.\t.\t.};

  $runner = Bio::EnsEMBL::VEP::Runner->new({
  %{$test_cfg->base_testing_cfg},
  input_data => $in,
  output_file => $test_cfg->{user_file}.'.out',
  gff => $test_cfg->{custom_gff},
  no_stats => 1,
  dont_skip => 1,
  quiet => 1,
  });

  $runner->run;
  open IN, $test_cfg->{user_file}.'.out';
  my @tmp_lines = <IN>;
  close IN;
  unlink($test_cfg->{user_file}.'.out');
  unlink($test_cfg->{user_file}.'.out_warnings.txt');
  is(scalar (grep {/SEQ_21/} @tmp_lines), 1, 'dont_skip variants which are not in GFF file');
}

done_testing();
