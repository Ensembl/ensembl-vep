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

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript');



## DATABASE TESTS
#################

SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 83 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_vepiens') if $can_use_db;

  # need to get a config object for further tests
  use_ok('Bio::EnsEMBL::VEP::Config');

  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    everything => 1,
    xref_refseq => 1,
  });
  ok($cfg, 'get new config object');
  
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $cfg
  });

  ok($as, 'new is defined');

  is(ref($as), 'Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript', 'check class');

  my $features;


  ## METHODS
  ##########

  is_deeply($as->valid_chromosomes, [21, 'LRG_485'], 'valid_chromosomes');

  is_deeply(
    $as->chr_lengths,
    {
      21 => 46709983,
      'LRG_485' => 9429
    },
    'chr lengths'
  );

  is_deeply(
    $as->info,
    {
      'genebuild' => '2014-07',
      'gencode' => 'GENCODE 24',
      'assembly' => 'GRCh38.p5',
      'polyphen' => '2.2.2',
      'sift' => 'sift5.2.2',
    },
    'info'
  );

  ok($as->check_sift_polyphen, 'check_sift_polyphen - not specified');

  $as->{sift} = 1;
  ok($as->check_sift_polyphen, 'check_sift_polyphen - SIFT');
  $as->{sift} = 0;

  $as->{polyphen} = 1;
  ok($as->check_sift_polyphen, 'check_sift_polyphen - PolyPhen');
  $as->{polyphen} = 0;

  my $ta = $multi->get_DBAdaptor('core')->get_TranscriptAdaptor;
  is(ref($ta), 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor', 'get transcript adaptor');
  my $tr = $ta->fetch_by_stable_id('ENST00000307301');
  is(ref($tr), 'Bio::EnsEMBL::Transcript', 'get transcript');

  ok($as->prefetch_gene_ids($tr), 'prefetch_transcript_ids');
  is($tr->{_gene_symbol}, 'MRPL39', 'prefetch_transcript_ids - gene symbol');
  is($tr->{_gene_symbol_source}, 'HGNC', 'prefetch_transcript_ids - gene symbol source');
  is($tr->{_gene_hgnc_id}, 'HGNC:14027', 'prefetch_transcript_ids - gene HGNC ID');

  ok($as->prefetch_transcript_ids($tr), 'prefetch_transcript_ids');
  is($tr->{_ccds}, 'CCDS33522.1', 'prefetch_transcript_ids - CCDS');
  is($tr->{_refseq}, 'NM_080794.3', 'prefetch_transcript_ids - RefSeq');

  ok($as->prefetch_translation_ids($tr, $tr->translation), 'prefetch_translation_ids');
  is($tr->{_swissprot}, 'Q9NYK5', 'prefetch_translation_ids - SWISSPROT');
  is($tr->{_uniparc}, 'UPI00001AEAC0', 'prefetch_translation_ids - UniParc');
  is($tr->{_protein}, 'ENSP00000305682', 'prefetch_translation_ids - Ensembl protein');

  ok($as->prefetch_translation_data($tr, $tr->translation), 'prefetch_translation_data');

  my $vep_cache = $tr->{_variation_effect_feature_cache};

  is(
    $vep_cache->{peptide},
    'MEALAMGSRALRLWLVAPGGGIKWRFIATSSASQLSPTELTEMRNDLFNKEKARQLSLTPRTEKIEVKHVGKTDPGTVFVMNKNISTPYSCAMHLSEWYCRKSILALVDGQPWDMYKPLTKSCEIKFLTFKDCDPGEVNKAYWRSCAMMMGCVIERAFKDEYMVNLVRAPEVPVISGAFCYDVVLDSKLDEWMPTKENLRSFTKDAHALIYKDLPFETLEVEAKVALEIFQHSKYKVDFIEEKASQNPERIVKLHRIGDFIDVSEGPLIPRTSICFQYEVSAVHNLQPTQPSLIRRFQGVSLPVHLRAHFTIWDKLLERSRKMTPFPILLLFTTQSFFTTSPESYLLHGTVSE',
    'prefetch_translation_data - peptide'
  );

  is_deeply(
    [sort {$a->{start} <=> $b->{start}} @{$vep_cache->{protein_features}}],
    [
      bless( {
        'end' => 324,
        'analysis' => bless( {
          '_display_label' => 'hmmpanther'
        }, 'Bio::EnsEMBL::Analysis' ),
        'hseqname' => 'PTHR11451',
        'start' => 38
      }, 'Bio::EnsEMBL::ProteinFeature' ),
      bless( {
        'end' => 132,
        'analysis' => bless( {
          '_display_label' => 'Gene3D'
        }, 'Bio::EnsEMBL::Analysis' ),
        'hseqname' => '3.10.20.30',
        'start' => 85
      }, 'Bio::EnsEMBL::ProteinFeature' ),
      bless( {
        'end' => 130,
        'analysis' => bless( {
          '_display_label' => 'Superfamily domains'
        }, 'Bio::EnsEMBL::Analysis' ),
        'hseqname' => 'SSF81271',
        'start' => 86
      }, 'Bio::EnsEMBL::ProteinFeature' ),
      bless( {
        'end' => 313,
        'analysis' => bless( {
          '_display_label' => 'Superfamily domains'
        }, 'Bio::EnsEMBL::Analysis' ),
        'hseqname' => 'SSF55186',
        'start' => 136
      }, 'Bio::EnsEMBL::ProteinFeature' ),
      bless( {
        'end' => 225,
        'analysis' => bless( {
          '_display_label' => 'Gene3D'
        }, 'Bio::EnsEMBL::Analysis' ),
        'hseqname' => '1tkeA02',
        'start' => 138
      }, 'Bio::EnsEMBL::ProteinFeature' ),
      bless( {
        'end' => 267,
        'analysis' => bless( {
          '_display_label' => 'Gene3D'
        }, 'Bio::EnsEMBL::Analysis' ),
        'hseqname' => '1tkeA03',
        'start' => 226
      }, 'Bio::EnsEMBL::ProteinFeature' ),
      bless( {
        'end' => 341,
        'analysis' => bless( {
          '_display_label' => 'Low complexity (Seg)'
        }, 'Bio::EnsEMBL::Analysis' ),
        'hseqname' => 'seg',
        'start' => 329
      }, 'Bio::EnsEMBL::ProteinFeature' ),
    ],
    'prefetch_translation_data - protein features'
  );

  my $tr2 = $ta->fetch_by_stable_id('ENST00000352957');
  $as->{$_} = 1 for qw(sift polyphen);
  $as->prefetch_translation_data($tr2, $tr2->translation);
  $vep_cache = $tr2->{_variation_effect_feature_cache};

  is(
    ref($vep_cache->{protein_function_predictions}->{sift}),
    'Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix',
    'prefetch_translation_data - sift matrix'
  );

  is(
    ref($vep_cache->{protein_function_predictions}->{polyphen_humvar}),
    'Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix',
    'prefetch_translation_data - polyphen_humvar matrix'
  );

  is(
    ref($vep_cache->{protein_function_predictions}->{polyphen_humdiv}),
    'Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix',
    'prefetch_translation_data - polyphen_humdiv matrix'
  );

  $vep_cache = $tr->{_variation_effect_feature_cache};

  ok($as->prefetch_transcript_data($tr), 'prefetch_transcript_data');
  is(ref($vep_cache->{introns}->[0]), 'Bio::EnsEMBL::Intron', 'prefetch_transcript_data - introns');
  is(ref($vep_cache->{sorted_exons}->[0]), 'Bio::EnsEMBL::Exon', 'prefetch_transcript_data - sorted_exons 1');
  is($vep_cache->{sorted_exons}->[0]->stable_id, 'ENSE00003528074', 'prefetch_transcript_data - sorted_exons 2');
  is(ref($vep_cache->{mapper}), 'Bio::EnsEMBL::TranscriptMapper', 'prefetch_transcript_data - mapper');
  is($vep_cache->{codon_table}, 1, 'prefetch_transcript_data - codon_table');
  
  is(
    $vep_cache->{five_prime_utr}->seq,
    'AGGGCGGAGAAGGACTTGCGCGCGACGGTTCTCACCGCTGCT',
    'prefetch_transcript_data - five_prime_utr'
  );
  
  is(
    $vep_cache->{three_prime_utr}->seq,
    'CTGAAGATCAAAGTAAAGCAACAGAGGAATGTACATCTACCTAATAACTTTCTAAAATTTAAATATGTATAATAAAATAAATGTTTTAAATATAA',
    'prefetch_transcript_data - three_prime_utr'
  );

  my $tr_same_utr = $ta->fetch_by_stable_id('ENST00000352957');
  ok($as->prefetch_transcript_data($tr_same_utr), 'same utr - prefetch_transcript_data');
  ok($vep_cache->{five_prime_utr} eq $tr_same_utr->{_variation_effect_feature_cache}->{five_prime_utr}, 'same utr - same ref');
  
  is(
    $vep_cache->{translateable_seq},
    'ATGGAGGCGCTGGCCATGGGTTCCCGGGCGCTGCGGCTCTGGCTGGTCGCACCCGGTGGCGGGATCAAATGGAGATTTATAGCAACATCGTCAGCTTCTCAGCTGTCACCGACAGAATTGACAGAAATGCGGAATGATCTCTTTAATAAAGAGAAAGCCAGGCAGTTATCATTAACTCCCCGAACTGAGAAGATAGAAGTTAAGCATGTTGGGAAAACTGACCCCGGTACTGTCTTCGTGATGAATAAAAACATTTCAACTCCCTACAGTTGTGCCATGCATTTAAGCGAGTGGTATTGCAGGAAGTCCATTCTGGCTCTGGTGGATGGACAGCCTTGGGACATGTATAAGCCTTTAACAAAGTCCTGTGAAATTAAATTTCTTACTTTCAAAGATTGTGATCCAGGAGAAGTGAATAAGGCATATTGGCGTTCCTGTGCTATGATGATGGGCTGTGTGATAGAGAGGGCATTCAAAGATGAATATATGGTCAATTTGGTCAGAGCTCCAGAAGTTCCTGTAATTTCTGGTGCCTTCTGTTATGACGTAGTTTTGGATAGCAAACTTGATGAGTGGATGCCAACAAAAGAGAACTTACGTTCCTTCACAAAAGATGCTCATGCTTTAATTTATAAAGATCTTCCATTTGAAACTCTGGAAGTTGAAGCAAAAGTGGCATTGGAAATATTTCAACACAGCAAGTACAAAGTAGATTTCATAGAAGAGAAGGCATCTCAGAACCCTGAGAGAATAGTCAAGCTACACAGAATAGGTGACTTCATTGATGTGAGTGAGGGCCCTCTTATTCCAAGAACAAGTATTTGTTTCCAGTATGAAGTATCAGCAGTTCACAATCTTCAACCCACCCAGCCAAGTCTCATACGAAGATTCCAGGGCGTGTCTTTACCTGTTCACTTAAGAGCACATTTTACAATATGGGATAAGCTATTGGAAAGATCTCGGAAAATGACTCCATTTCCCATTCTCCTTCTATTTACTACGCAGTCATTCTTCACTACCTCGCCTGAGTCGTACCTCCTCCATGGAACAGTCTCAGAGTAA',
    'prefetch_transcript_data - translateable_seq'
  );
  

  $features = $as->get_features_by_regions_uncached([[21, 511]]);
  is(ref($features), 'ARRAY', 'get_features_by_regions_uncached ref 1');
  is(scalar @$features, 6, 'get_features_by_regions_uncached count');
  is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'get_features_by_regions_uncached ref 2');
  is($features->[0]->stable_id, 'ENST00000456917', 'get_features_by_regions_uncached stable_id');

  # now we should be able to retrieve the same from memory
  $features = $as->get_features_by_regions_cached([[21, 511]]);
  is(ref($features), 'ARRAY', 'get_features_by_regions_cached ref 1');
  is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'get_features_by_regions_cached ref 2');
  is($features->[0]->stable_id, 'ENST00000456917', 'get_features_by_regions_cached stable_id');

  # test merged flag setting transcript source
  $as->{merged} = 1;
  $features = $as->get_features_by_regions_uncached([[21, 511]]);
  is($features->[0]->{_source_cache}, 'Ensembl', 'get_features_by_regions_uncached - merged sets _source_cache');
  $as->{merged} = undef;

  $as->clean_cache();
  is_deeply($as->cache, {}, 'clean_cache');



  ## TEST REFSEQ
  ##############

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    refseq => 1,
  });
  ok($cfg, 'refseq - get new config object');
  
  $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $cfg
  });

  $ta = $multi->get_DBAdaptor('otherfeatures')->get_TranscriptAdaptor;
  is(ref($ta), 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor', 'refseq - get transcript adaptor');
  $tr = $ta->fetch_by_stable_id('NM_201413.2');
  is(ref($tr), 'Bio::EnsEMBL::Transcript', 'refseq - get transcript');

  $as->prefetch_gene_ids($tr);
  $as->prefetch_transcript_ids($tr);

  is($tr->{_gene_symbol}, 'APP', 'refseq - prefetch_transcript_ids - gene symbol');

  $features = $as->get_features_by_regions_uncached([[21, 511]]);
  is(ref($features), 'ARRAY', 'refseq - get_features_by_regions_uncached ref 1');
  is(scalar @$features, 8, 'refseq - get_features_by_regions_uncached count');
  is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'refseq - get_features_by_regions_uncached ref 2');
  is($features->[0]->stable_id, 'NR_001458.3', 'refseq - get_features_by_regions_uncached stable_id');

  # now we should be able to retrieve the same from memory
  $features = $as->get_features_by_regions_cached([[21, 511]]);
  is(ref($features), 'ARRAY', 'refseq - get_features_by_regions_cached ref 1');
  is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'refseq - get_features_by_regions_cached ref 2');
  is($features->[0]->stable_id, 'NR_001458.3', 'refseq - get_features_by_regions_cached stable_id');

  $as->clean_cache();



  ## TESTS WITH AN INPUT BUFFER
  #############################

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    everything => 1,
    xref_refseq => 1,
  });
  
  $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $cfg
  });


  use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
  my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
  ok($p, 'get parser object');

  use_ok('Bio::EnsEMBL::VEP::InputBuffer');
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  is(ref($ib), 'Bio::EnsEMBL::VEP::InputBuffer', 'check class');

  is(ref($ib->next()), 'ARRAY', 'check buffer next');

  is_deeply(
    $as->get_all_regions_by_InputBuffer($ib),
    [
      [
        '21',
        511
      ],
      [
        '21',
        512
      ],
      [
        '21',
        513
      ],
      [
        '21',
        514
      ],
      [
        '21',
        515
      ],
      [
        '21',
        517
      ],
      [
        '21',
        518
      ],
      [
        '21',
        519
      ]
    ],
    'get_all_regions_by_InputBuffer'
  );

  $features = $as->get_all_features_by_InputBuffer($ib);
  is(ref($features), 'ARRAY', 'get_all_features_by_InputBuffer ref 1');
  is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'get_all_features_by_InputBuffer ref 2');
  is(ref($features->[-1]), 'Bio::EnsEMBL::Transcript', 'get_all_features_by_InputBuffer ref 3');
  is($features->[0]->stable_id, 'ENST00000567517', 'get_all_features_by_InputBuffer stable_id');
  is(scalar @$features, 44, 'get_all_features_by_InputBuffer count');

  # do it again to get them from memory
  $features = $as->get_all_features_by_InputBuffer($ib);
  is($features->[0]->stable_id, 'ENST00000567517', 'get_all_features_by_InputBuffer again');

  $ib->next();
  is_deeply($as->get_all_features_by_InputBuffer($ib), [], 'get_all_features_by_InputBuffer on empty buffer');

  # test upstream
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->create_input_file([qw(21 13002936 . C A)])});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  is_deeply(
    $as->get_all_regions_by_InputBuffer($ib),
    [
      [
        '21',
        259
      ],
      [
        '21',
        260
      ],
    ],
    'get_all_regions_by_InputBuffer - include upstream'
  );


  # reset
  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  $as->annotate_InputBuffer($ib);
  my $vf = $ib->buffer->[0];
  my $tvs = $vf->get_all_TranscriptVariations;

  is(scalar @$tvs, 3, 'annotate_InputBuffer - get_all_TranscriptVariations count');

  $vf->_finish_annotation;
  is($vf->display_consequence, 'missense_variant', 'annotate_InputBuffer - display_consequence');


  ## transcript filter
  $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $cfg,
    filter => 'stable_id ne ENST00000352957',
  });

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}, valid_chromosomes => [21]});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  $as->annotate_InputBuffer($ib);
  $vf = $ib->buffer->[0];
  $vf->_finish_annotation;
  is(scalar (grep {$_->transcript->stable_id eq 'ENST00000352957'} @{$vf->get_all_TranscriptVariations}), 0, 'with filter - filtered transcript absent');
  is($vf->display_consequence, '3_prime_UTR_variant', 'with filter - display_consequence');



  ## TEST REFSEQ
  ##############

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'homo_vepiens',
    refseq => 1,
  });
  
  $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $cfg
  });

  is_deeply(
    $as->info,
    {
      'refseq' => '2015-10-26 11:20:17 - GCF_000001405.31_GRCh38.p5_genomic.gff',
      'genebuild' => '2014-07',
      'gencode' => 'GENCODE 24',
      'assembly' => 'GRCh38.p5'
    },
    'refseq - info'
  );

  $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $test_cfg->{test_vcf}});
  $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();

  $features = $as->get_all_features_by_InputBuffer($ib);
  is(ref($features), 'ARRAY', 'refseq - get_all_features_by_InputBuffer ref 1');
  is(ref($features->[0]), 'Bio::EnsEMBL::Transcript', 'refseq - get_all_features_by_InputBuffer ref 2');
  is(ref($features->[-1]), 'Bio::EnsEMBL::Transcript', 'refseq - get_all_features_by_InputBuffer ref 3');
  is($features->[0]->stable_id, 'NR_024092.1', 'refseq - get_all_features_by_InputBuffer stable_id');
  is(scalar @$features, 35, 'refseq - get_all_features_by_InputBuffer count');

  $as->clean_cache();
  $as->{all_refseq} = 1;
  $features = $as->get_all_features_by_InputBuffer($ib);
  is(scalar @$features, 51, 'refseq - get_all_features_by_InputBuffer count - all_refseq');
  $as->{all_refseq} = 0;



  ## TEST ANOTHER SPECIES
  #######################

  $multi = Bio::EnsEMBL::Test::MultiTestDB->new('mus_muscuvep');

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'mus_muscuvep',
    everything => 1,
    quiet => 1,
  });
  ok($cfg, 'get new config object');

  ok($cfg->param('sift'), 'other species - sift on before setup with --everything');
  ok($cfg->param('polyphen'), 'other species - polyphen on before setup with --everything');
  
  # instantiating this should switch off --sift and --polyphen that are passively switched on by --everything
  $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $cfg
  });

  ok(!$cfg->param('sift'), 'other species - sift off after setup with --everything');
  ok(!$cfg->param('polyphen'), 'other species - polyphen off after setup with --everything');

  $cfg = Bio::EnsEMBL::VEP::Config->new({
    %$db_cfg,
    database => 1,
    offline => 0,
    species => 'mus_muscuvep',
    sift => 'b',
    quiet => 1,
  });

  throws_ok {
    $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $cfg
  })} qr/SIFT not available/, 'other species - throws on PolyPhen request';

  1;
}

# done
done_testing();
