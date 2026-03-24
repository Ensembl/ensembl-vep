use Bio::EnsEMBL::VEP::Runner;
use Bio::EnsEMBL::VEP::OutputFactory::JSON;
use Data::Dumper;

use FindBin qw($Bin);
use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

sub get_annotated_buffer {
  my $tmp_cfg = shift;

  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    %$cfg_hash,
    dir => $test_cfg->{cache_root_dir},
    %$tmp_cfg,
  });

  $runner->init;

  my $ib = $runner->get_InputBuffer;
  $ib->next();
  $_->annotate_InputBuffer($ib) for @{$runner->get_all_AnnotationSources};
  $ib->finish_annotation();

  return $ib;
}

my $ib = get_annotated_buffer({
    input_file => $test_cfg->{test_vcf},
    everything => 1,
    dir => $test_cfg->{cache_root_dir},
  });


$of = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $ib->config});
@lines = @{$of->get_all_lines_by_InputBuffer($ib)};

print Dumper \@lines, "\n";
exit 0;
############################
# use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

# our @SORTED_OVERLAP_CONSEQUENCES = sort {$a->tier <=> $b->tier} values %OVERLAP_CONSEQUENCES;

# foreach (@SORTED_OVERLAP_CONSEQUENCES) {
#     print $_->{SO_term}, " ", $_->{tier}, "\n";
# }

use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Digest::MD5 qw(md5_hex);

my $registry = 'Bio::EnsEMBL::Registry';
# #$registry->load_registry_from_db(
# #    -host => 'mysql-ens-var-prod-4',
# #    -user => 'ensro',
# #    -port => 4694,
# #    -pass => '',
# #    -db_version   => 116,
# #    
# #);
$registry->load_all("/hps/nobackup/flicek/ensembl/variation/snhossain/release/116/human/protein_function_dbnsfp_test/ensembl.registry");
$registry->add_alias('homo_sapiens', 'human');

# my $transcript_adaptor = $registry->get_adaptor('homo_sapiens', 'Core', 'Transcript');
# my $pfpm_adaptor = $registry->get_adaptor( 'human', 'Variation', 'ProteinFunctionPredictionMatrix' );
# my $transcript = $transcript_adaptor->fetch_by_stable_id("ENST00000311351");
# print $transcript->stable_id, "\n";
# my $translation = $transcript->translation;
# my $md5 = md5_hex($translation->seq);
# $md5="6a9416abf1715498230a15394688bcce";
# $pfpm = $pfpm_adaptor->fetch_by_analysis_translation_md5( "dbnsfp_esm1b", $md5 );
# if (defined $pfpm) {
#     my $matrix = $pfpm->deserialize;
#     # print Dumper($matrix), "\n" if defined $pfpm;
#     print Dumper $matrix->{86}->{A}, "\n";
#     # foreach my $pos (keys %{$matrix}) {
#     #     foreach my $aa (%{ $matrix->{$pos} }) {
#     #             my ($pred, $score) = @{$matrix->{$pos}->{$aa}};
#     #             print($pos, ":", $aa, "-", $pred, " ", $score, "\n");
#     #     }
#     # }
# }
# exit(0)



use Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation;
use Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation;
use Bio::EnsEMBL::Variation::DBSQL::AttributeAdaptor;
# use Bio::EnsEMBL::Registry;
# use Data::Dumper;
# use Digest::MD5 qw(md5_hex);
# use List::MoreUtils qw(firstidx);

# my $registry = 'Bio::EnsEMBL::Registry';
# $registry->load_all("/hps/nobackup/flicek/ensembl/variation/snhossain/release/116/human/protein_function_dbnsfp_test/ensembl.registry");

# my $transcript_adaptor = $registry->get_adaptor("homo_sapiens", "core", "transcript");
# # my $tr = $transcript_adaptor->fetch_by_stable_id("ENST00000311351");
# my $tr = $transcript_adaptor->fetch_by_stable_id("rna-NM_014462.3");
# my $tl = $tr->translation;
# print $tl->stable_id, "\n";
# print md5_hex($tl->seq), "\n";
# exit 0;

# my $dbNSFP = Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation->new(
#     -registry_file => '/hps/nobackup/flicek/ensembl/variation/snhossain/release/116/human/protein_function_dbnsfp_test/ensembl.registry',
#     -species => 'homo_sapiens',
#     -working_dir => '/hps/software/users/ensembl/variation/snhossain',
#     -annotation_file =>  '/nfs/production/flicek/ensembl/variation/data/dbNSFP/5.2a/dbNSFP5.2a_grch38.gz',
#     -assembly => 'GRCh38',
#     -annotation_file_version => '5.2a',
#   );
# # $dbNSFP->run('756a8a511d1f4557d4e577fa7d55df1a');
# $dbNSFP->run('6a9416abf1715498230a15394688bcce', {'6a9416abf1715498230a15394688bcce' => 'ENSP00000310596'});
# exit(0)

my $cadd = Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation->new(
    -registry_file => '/hps/nobackup/flicek/ensembl/variation/snhossain/release/116/human/protein_function_dbnsfp_test/ensembl.registry',
    -species => 'homo_sapiens',
    -working_dir => '/hps/software/users/ensembl/variation/snhossain',
    -annotation_file =>  '/nfs/production/flicek/ensembl/variation/data/CADD/v1.7/grch38/CADD_GRCh38_1.7_whole_genome_SNVs.tsv.gz',
    -assembly => 'GRCh38',
    -annotation_file_version => 'v1.7',
  );
$cadd->run('6a9416abf1715498230a15394688bcce', {'6a9416abf1715498230a15394688bcce' => 'ENSP00000310596'});
exit(0)

# sub get_gvf_line {
#   my ($line, $id_count) = @_;
#   my $gvf_line = {};
#   my @header_names = qw/seq_id source type start end score strand phase/;
#   my @header_values = split(/\t/, $line);
#   my $attrib = pop @header_values; 

#   for my $i (0 .. $#header_names) { 
#     $gvf_line->{$header_names[$i]} = $header_values[$i]; 
#   }

#   my @attributes = split(';', $attrib);
#   foreach my $attribute (@attributes) {
#     my ($key, $value) = split('=', $attribute);
#     if ($value) {
#       $gvf_line->{attributes}->{$key} = $value;
#     }
#   }

#   $gvf_line->{attributes}->{ID} = $id_count;
#   $line = join("\t", map {$gvf_line->{$_}} (
#     'seq_id', 
#     'source', 
#     'type', 
#     'start', 
#     'end', 
#     'score', 
#     'strand', 
#     'phase'));
#   my $attributes = join(";", map{"$_=$gvf_line->{attributes}->{$_}"} keys %{$gvf_line->{attributes}});
#   return "$line\t$attributes";
# }

# my $dir = "/hps/nobackup/flicek/ensembl/variation/snhossain/release/114/ovis_aries/release_dump/vertebrates/variation/gvf/ovis_aries_rambouillet";
# my $dump_type = "ovis_aries_rambouillet_incl_consequences";
# my $fh_join = FileHandle->new("$dir/$file_name.gvf", 'a');
# my $input_ids = [13,14,15,16,17,18,19,2,'22_21',23,27,'28_25',3,'39_20',4,5,6,7,8,9];
# foreach my $file_id (@$input_ids) {
#     if (-e "$dir/$dump_type-$file_id.gvf.gz") {
#       `gunzip $dir/$dump_type-$file_id.gvf.gz`; 
#     }
#     my $fh = FileHandle->new("$dir/$dump_type-$file_id.gvf", 'r');
#     print "$dir/$dump_type-$file_id.gvf", "\n";
#     while (<$fh>) {
#       chomp;
#       my $line = $_;
#       next if ($line =~ m/^#/);
#       my $gvf_line = get_gvf_line($line, $id_count);
#       print $fh_join $gvf_line, "\n";   
#       $id_count++;
#     }
#     $fh->close();

#     # `rm $dir/$dump_type-$file_id.gvf`;
#   }
# $fh_join->close();