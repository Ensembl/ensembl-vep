=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::QCDump;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::BaseVEP);

use File::Path qw(make_path rmtree);
use Digest::MD5 qw(md5_hex);
use JSON;
use Scalar::Util qw(looks_like_number);
use FileHandle;
use Data::Dumper;
$Data::Dumper::Indent = 1;

use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::CacheDir;
use Bio::EnsEMBL::VEP::Runner;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

sub param_defaults {
  return {
    'variation'  => 0,
    'regulation' => 0,
    'dir_suffix' => '',
    'convert'    => 0,
  };
}

sub run {
  my $self = shift;

  my @report_files;

  push @report_files, @{$self->qc()};
  push @report_files, @{$self->qc('_tabixconverted')} if $self->param('convert') && $self->param('variation');

  if(@report_files) {
    die(
      "ERROR: Encountered some errors in test set checks, see the following files for details:\n".
      join("\n", @report_files)
    );
  }

  return 1;
}

sub qc {
  my ($self, $mod) = @_;

  my $type      = $self->param('type');
  my $has_var   = $self->param('variation');
  my $has_reg   = $self->param('regulation');
  my $converted = $mod && $mod =~ /tabix/;
  my $species   = $self->required_param('species');
  my $assembly  = $self->required_param('assembly');

  my $dump_dir = $self->dump_dir();
  my $tar_file = $self->get_tar_file_name(
    $dump_dir,
    $species.($type eq 'core' ? '' : "_$type"),
    $assembly,
    $mod
  );

  die("ERROR: Tar file $tar_file not found\n") unless -e $tar_file;

  my $qc_dir = $dump_dir.'/qc/'.md5_hex($self->input_id.($converted || 0));

  unless(-d $qc_dir) {
    make_path($qc_dir) or die "Could not make directory $qc_dir";
    $self->run_system_command(sprintf('tar -C %s -xzf %s', $qc_dir, $tar_file));
  }

  # check dir exists
  my $method_name = ($type eq 'core' ? '' : $type.'_').'species_suffix';
  my $extracted_dir = $qc_dir.'/'.$self->$method_name;
  die("ERROR: Expected to find $extracted_dir\n") unless -d $extracted_dir;

  # check info.txt exists
  die("ERROR: Expected to find $extracted_dir/info.txt") unless -e "$extracted_dir/info.txt";

  # create objects
  my $config_obj = Bio::EnsEMBL::VEP::Config->new({
    dir => $extracted_dir,
    offline => 1,
    species => $species,
    assembly => $assembly,
    check_existing => 1,
    regulatory => 1,
  });
  my $cache_dir_obj = Bio::EnsEMBL::VEP::CacheDir->new({dir => $extracted_dir, config => $config_obj});

  # these subs check different aspects of the cache
  $self->check_info($cache_dir_obj, $converted);

  $self->check_annotation_sources($cache_dir_obj, $converted);

  $self->check_dirs($cache_dir_obj, $converted);

  my @report_files = ();
  my $report;

# Don't run for release/94 this keeps causing errors ENSVAR-1172
# $report = $self->human_frequency_checks($qc_dir, $mod) if $has_var && $species eq 'homo_sapiens';
# push @report_files, $report if $report;

#  $report = $self->run_test_set($qc_dir, $mod) if $has_var && $type eq 'core';
#  push @report_files, $report if $report;

  # clean up
  rmtree($qc_dir);

  return \@report_files;
}

sub check_info {
  my ($self, $cache_dir_obj, $converted) = @_;

  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');

  my $info = $cache_dir_obj->info;
  die("ERROR: No keys in info\n") unless keys %$info;

  die("ERROR: species info key does not exist\n") unless $info->{species};
  die("ERROR: species info key does not match job species name\n")
    unless $info->{species} eq $species;

  die("ERROR: assembly info key does not exist\n") unless $info->{assembly};
  die("ERROR: assembly info key does not match job assembly name\n")
    unless $info->{assembly} eq $assembly;

  # var_type should be tabix if converted
  if($converted) {
    die("ERROR: var_type info key is not set to tabix\n")
      unless $info->{var_type} && $info->{var_type} eq 'tabix';
  }

  # check version_data
  die("ERROR: version_data key not found\n") unless $info->{version_data} && ref($info->{version_data}) eq 'HASH';

  # check variation_cols defined
  if($self->param('variation')) {
    die("ERROR: variation_cols info key not found\n") unless $info->{variation_cols};

    die("ERROR: variation_cols info value looks wrong: ".$info->{variation_cols}."\n")
      unless ref($info->{variation_cols}) eq 'ARRAY';

    if($converted) {
      die("ERROR: variation_cols value does not include chr\n")
        unless $info->{variation_cols}->[0] eq 'chr';
    }
  }
}

sub check_annotation_sources {
  my ($self, $cache_dir_obj, $converted) = @_;

  my $annotation_sources = $cache_dir_obj->get_all_AnnotationSources;

  die("ERROR: No annotation sources retrieved\n") unless ref($annotation_sources) eq 'ARRAY' && scalar @$annotation_sources;

  # expect core always
  die("ERROR: No Transcript annotation source found\n")
    unless grep {ref($_) eq 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript'} @$annotation_sources;

  if($self->param('variation')) {
    my $t = 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation'.($converted ? 'Tabix' : '');
    die("ERROR: No Variation annotation source found\n")
      unless grep {ref($_) eq $t} @$annotation_sources;
  }
  if($self->param('regulation')) {
    die("ERROR: No Regulation annotation source found\n")
      unless grep {ref($_) eq 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat'} @$annotation_sources;
  }
}

sub check_dirs {
  my ($self, $cache_dir_obj, $converted) = @_;

  my $dir = $cache_dir_obj->dir;
  opendir DIR, $dir or die $!;
  my %dirs_in_cache;
  $dirs_in_cache{$_} = 1 for grep {!/^\./ && -d "$dir/$_" } readdir(DIR);
  closedir DIR;

  die("ERROR: No chromosome dirs found\n") unless scalar keys %dirs_in_cache;

  # connect to DB to find all the slices we should have dirs for
  Bio::EnsEMBL::Registry->clear();
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -group   => 'core',
    -species => $self->param('species'),
    -port    => $self->param('port'),
    -host    => $self->param('host'),
    -user    => $self->param('user'),
    -pass    => $self->param('pass'),
    -dbname  => $self->param('dbname'),
    -MULTISPECIES_DB => $self->param('is_multispecies'),
    -species_id => $self->param('species_id'),
  );

  # get slices
  my $sa = $dba->get_SliceAdaptor;
  my @slices = @{$sa->fetch_all('toplevel')};
  push @slices, map {$_->alternate_slice} map {@{$_->get_all_AssemblyExceptionFeatures}} @slices;
  push @slices, @{$sa->fetch_all('lrg', undef, 1, undef, 1)} if $self->param('lrg');

  # check for dirs that don't appear in this list of slices
  delete($dirs_in_cache{$_->seq_region_name}) for @slices;
  die("ERROR: Found ".(scalar keys %dirs_in_cache)." dirs in cache that have no corresponding slice\n") if scalar keys %dirs_in_cache;

  # now filter out slices with no transcripts
  my $ta = $dba->get_TranscriptAdaptor();
  @slices = grep {$ta->count_all_by_Slice($_)} @slices;

  my $has_var = $self->param('variation') ? 1 : 0;
  my $has_reg = $self->param('regulation') ? 1 : 0;

  foreach my $slice(@slices) {
    my $chr = $slice->seq_region_name;
    die("ERROR: Dir $dir/$chr missing\n") unless -d "$dir/$chr";

    opendir SUB, "$dir/$chr" or die $!;
    my @files = grep {!/^\./} readdir(SUB);
    closedir SUB;

    die("ERROR: $dir/$chr is empty\n") unless @files;

    foreach my $file(@files) {
      die("ERROR: File name $file does not look right\n")
        unless $file =~ /\d+\-\d+(\_(var|reg))?\.gz/
        or $file eq 'all_vars.gz'
        or $file eq 'all_vars.gz.tbi'
        or $file eq 'all_vars.gz.csi';

      if($file =~ /var/) {
        die("ERROR: Cache should not contain var files, found $file in $dir/$chr\n") unless $has_var;

        if($converted) {
          die("ERROR: Converted cache should not contain per-MB files, found $file in $dir/$chr\n") if $file =~ /\_var\.gz/;

          if($file eq 'all_vars.gz') {
            die("ERROR: Index file for $dir/$chr/$file missing") unless -e "$dir/$chr/all_vars.gz.tbi" or -e "$dir/$chr/all_vars.gz.csi";
          }
          elsif($file =~ /^all_vars\.gz\.(cs|tb)i$/) {
            die("ERROR: Data file for $dir/$chr/$file missing") unless -e "$dir/$chr/all_vars.gz";
          }
        }
        else {
          die("ERROR: Basic cache should not contain $file file, found $file in $dir/$chr\n") if $file =~ /^all_vars/;
        }
      }
      elsif($file =~ /reg/) {
        die("ERROR: Cache should not contain regulation files, found $file in $dir/$chr\n") unless $has_reg;
      }
    }
  }
}

sub run_test_set {
  my $self = shift;
  my $qc_dir = shift;
  my $mod = shift;

  # get a test file
  my $test_file = $self->generate_test_file();

  my $has_var = $self->param('variation') ? 1 : 0;
  my $has_reg = $self->param('regulation') ? 1 : 0;

  Bio::EnsEMBL::Registry->clear();
  # create a VEP Runner
  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    species => $self->param('species'),
    assembly => $self->param('assembly'),
    db_version => $self->param('ensembl_release'),
    cache_version => $self->param('eg_version') || $self->param('ensembl_release'),
    offline => 0,
    cache => 1,
    database => 0,
    host => $self->param('host'),
    user => $self->param('user'),
    port => $self->param('port'),
    password => $self->param('pass'),
    is_multispecies => $self->param('is_multispecies'),
    dir => $qc_dir,
    input_file => $test_file,
    format => 'ensembl',
    output_format => 'json',
    safe => 1,
    quiet => 1,
    no_stats => 1,
    regulatory => $has_reg,
    check_existing => $has_var,
    buffer_size => 10,
    # check_ref => 1,
    warning_file => $qc_dir.'/warnings.txt',
    $self->param('type') => 1,
  });

  my $json = JSON->new();

  # qc report file
  my $qc_report_file = sprintf(
    '%s/qc/%s_%s_%s%s_qc_report.txt',
    $self->dump_dir(),
    $self->param('species'),
    $self->param('type'),
    $self->param('assembly'),
    $mod || "",
  );

  open QC, ">$qc_report_file" or die $!;

  $| = 1;

  while(my $line = $runner->next_output_line) {
    my $data = $json->decode($line);

    print QC "ERROR: input field not found in JSON output\nQC dir: $qc_dir\n".(Dumper $data)."\n\n" unless $data->{input};
    my @input = split("\t", $data->{input});

    my $feature_type = $input[-1];
    my $expected_cons = join(",", sort split(",", $input[-2]));
    my $feature_id = $input[-3];

    # check consequence type
    print QC "ERROR: no data for $feature_type found in JSON output\nQC dir: $qc_dir\n".(Dumper $data)."\n\n" unless $data->{$feature_type.'_consequences'};

    if(
      my ($blob) = grep {
        # motif_feature stable_id from DB is the regfeat ID, annoyingly
        $feature_type eq 'motif_feature' ? 1 : $_->{$feature_type.'_id'} eq $feature_id
      } @{$data->{$feature_type.'_consequences'}}
    ) {
      print QC "ERROR: no consequence_terms field in blob\nQC dir: $qc_dir\n".(Dumper $data)."\n\n" unless $blob->{consequence_terms};
      my $got_cons = join(",", sort @{$blob->{consequence_terms}});
      print QC "ERROR: consequence_types don't match, expected: $expected_cons, got: $got_cons\nQC dir: $qc_dir\n".(Dumper $data)."\n\n" unless $expected_cons eq $got_cons;
    }
    else {
      print QC "ERROR: no data for $feature_id found in JSON output\nQC dir: $qc_dir\n".(Dumper $data)."\n\n";
    }

    # check colocated variants
    if($has_var) {
      my $expected_var_id = $input[-4];
      print QC "ERROR: no data for colocated_variants found in JSON output\nQC dir: $qc_dir\n".(Dumper $data)."\n\n" unless $data->{colocated_variants};

      print QC "ERROR: expected var id $expected_var_id not found\nQC dir: $qc_dir\n".(Dumper $data)."\n\n"
        unless grep {$_->{id} eq $expected_var_id} @{$data->{colocated_variants}};
    }
  }

  close QC;

  my $count = `wc -l $qc_report_file`;

  # report file empty, no failures
  if($count =~ /^0/) {
    unlink($qc_report_file);
    return undef;
  }
  else {
    return $qc_report_file;
  }
}

sub generate_test_file {
  my $self = shift;

  my $file = sprintf(
    '%s/qc/%s_%s_test_input.txt',
    $self->dump_dir(),
    $self->param('species'),
    $self->param('assembly')
  );

  # use a lock file in case two processes try writing to the same one
  my $lock = $file.'.lock';

  return $file if -e $file && !-e $lock;

  # wait 3hrs for other process to finish, it will be a long-running query in human
  my $sleep_count = 0;
  if(-e $lock) {
    while(-e $lock) {
      sleep 1;
      die("I've been waiting for $lock to be removed for $sleep_count seconds, something may have gone wrong\n") if ++$sleep_count > (60*60*3);
    }

    die("ERROR: Finished waiting but $file not found\n") unless -e $file;

    return $file;
  }

  open OUT, ">$lock" or die $!;
  print OUT "$$\n";
  close OUT;
  $self->{_lock_file} = $lock;

  my $core_dbname = $self->param('dbname');
  $core_dbname =~ s/otherfeatures/core/;
  my $var_dbname = $core_dbname;
  $var_dbname =~ s/core/variation/;

  my $dba = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
    -group   => 'variation',
    -species => $self->param('species'),
    -port    => $self->param('port'),
    -host    => $self->param('host'),
    -user    => $self->param('user'),
    -pass    => $self->param('pass'),
    -dbname  => $var_dbname,
    -MULTISPECIES_DB => $self->param('is_multispecies'),
    -species_id => $self->param('species_id'),
  );

  # get distinct seq_region_ids
  # some species have huge numbers of seq_regions, so we're going to limit the number we analyse
  # by picking only the 50 with the most variation_features on
  my $sth = $dba->dbc->prepare(qq{
    SELECT seq_region_id, COUNT(*)
    FROM variation_feature
    GROUP BY seq_region_id
    ORDER BY COUNT(*) DESC
    LIMIT 50
  });
  $sth->execute;
  my @sr_ids;
  while(my $arrayref = $sth->fetchrow_arrayref) {
    push @sr_ids, $arrayref->[0];
  }
  $sth->finish;

  my @tables = ('transcript', ($self->param('regulation') ? qw(regulatory_feature motif_feature) : ()));

  my @rows;

  foreach my $table(@tables) {
    $sth = $dba->dbc->prepare(qq{
      SELECT distinct(consequence_types)
      FROM $table\_variation
    });
    $sth->execute;

    # The idea of this query is to fetch one input variant to test
    # for every conseqeunce_type and seq_region combination.
    # It contains a few extra clauses to:
    # - exclude failed variants
    # - exclude variants with non-standard allele strings (e.g. COSMIC)
    # - exclude transcript_variation entries where the variation_feature is not on the same seq_region as the transcript
    my $sth2 = $dba->dbc->prepare(qq{
      SELECT s.name, vf.seq_region_start, vf.seq_region_end, vf.allele_string, vf.seq_region_strand,
      vf.variation_name, tv.feature_stable_id, tv.consequence_types
      FROM seq_region s, $table\_variation tv, $core_dbname\.transcript t, variation_feature vf
      LEFT JOIN failed_variation fv on vf.variation_id = fv.variation_id
      WHERE s.seq_region_id = vf.seq_region_id
      AND vf.variation_feature_id = tv.variation_feature_id
      AND tv.feature_stable_id = t.stable_id
      AND vf.seq_region_id = t.seq_region_id
      AND fv.variation_id IS NULL
      AND vf.allele_string RLIKE '^[ACGT]+/[ACGT]+\$'
      AND tv.consequence_types = ?
      AND vf.seq_region_id = ?
      LIMIT 1
    });

    while(my $consref = $sth->fetchrow_arrayref) {
      foreach my $sr(@sr_ids) {
        $sth2->execute($consref->[0], $sr);

        while(my $arrayref = $sth2->fetchrow_arrayref) {
          push @rows, [@$arrayref, $table];
        }
      }
    }

    $sth->finish;
    $sth2->finish;
  }

  # print sorted in chr order
  my $fh = FileHandle->new("> $file");
  print $fh join("\t", @$_)."\n" for
    sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]}
    @rows;

  $fh->close();

  unlink($lock);

  return $file;
}

sub human_frequency_checks {
  my $self = shift;
  my $qc_dir = shift;
  my $mod = shift;

  # get a test file
  my $test_file = $self->generate_human_frequency_test_file();

  # read in data from it
  open IN, $test_file;
  my $ref_data_hash = {};
  while(<IN>) {
    chomp;
    my @tmp = split(/\s+/, $_);
    $ref_data_hash->{$tmp[5]} = [map {s/^.+?\://g; $_} @tmp];
  }
  close IN;

  my $has_var = $self->param('variation') ? 1 : 0;
  my $has_reg = $self->param('regulation') ? 1 : 0;

  Bio::EnsEMBL::Registry->clear();
  # create a VEP Runner
  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    species => $self->param('species'),
    assembly => $self->param('assembly'),
    db_version => $self->param('ensembl_release'),
    cache_version => $self->param('eg_version') || $self->param('ensembl_release'),
    offline => 0,
    cache => 1,
    database => 0,
    host => $self->param('host'),
    user => $self->param('user'),
    port => $self->param('port'),
    password => $self->param('pass'),
    is_multispecies => $self->param('is_multispecies'),
    dir => $qc_dir,
    input_file => $test_file,
    format => 'ensembl',
    delimiter => " ",
    output_format => 'tab',
    safe => 1,
    quiet => 1,
    no_stats => 1,
    check_existing => $has_var,
    af_1kg => 1,
    af_esp => 1,
    af_gnomad => 1,
    fields => 'Uploaded_variation,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF',
    buffer_size => 1,
    pick => 1,
    failed => 1,
    # check_ref => 1,
    warning_file => $qc_dir.'/warnings.txt',
    $self->param('type') => 1,
  });

  my $json = JSON->new();

  # qc report file
  my $qc_report_file = sprintf(
    '%s/qc/%s_%s_%s%s_frequency_qc_report.txt',
    $self->dump_dir(),
    $self->param('species'),
    $self->param('type'),
    $self->param('assembly'),
    $mod || "",
  );

  open QC, ">$qc_report_file" or die $!;

  $| = 1;

  while(my $line = $runner->next_output_line()) {
    my @data = split("\t", $line);
    my $rs = shift @data;

    if(my $ref_data = $ref_data_hash->{$rs}) {

      # ref_data contains all the input bits too, shift them off
      while(@$ref_data > @data) {
        shift @$ref_data;
      }

      # now compare
      my $mismatches = 0;

      for my $i(0..$#data) {
        next unless looks_like_number($ref_data->[$i]);
        $mismatches++ if !looks_like_number($data[$i]) || sprintf("%.3g", $ref_data->[$i]) != sprintf("%.3g", $data[$i]);
      }

      if($mismatches) {
        print QC "ERROR: Mismatched frequencies in $rs (IN vs OUT):\n".join("\t", @$ref_data)."\n".join("\t", @data)."\n";
      }
    }
    else {
      print QC "ERROR: no ref data found for $rs\n";
    }
  }

  close QC;

  my $count = `wc -l $qc_report_file`;

  # report file empty, no failures
  if($count =~ /^0/) {
    unlink($qc_report_file);
    return undef;
  }
  else {
    return $qc_report_file;
  }
}

sub generate_human_frequency_test_file {
  my $self = shift;

  my $assembly = $self->param('assembly');
  my $file = sprintf(
    '%s/qc/%s_%s_human_frequency_test_input.txt',
    $self->dump_dir(),
    $self->param('species'),
    $assembly
  );

  return $file if -e $file;

  open OUT, ">$file" or die $!;

## NB HARDCODED TEST VARIANTS
## These variants were randomly selected using a process analagous to the following:
## zcat ~/.vep/homo_sapiens/91_GRCh38/X/all_vars.gz | grep -v "," | grep -v "\/.*\/" | grep rs | perl select_vars.pl | recut  1,5,5,7,8,2,14- | sed "s/\t\.\t/\t\+\t/" | rand_lines.pl 0.01 | head | tr "\t" " "

## select_vars.pl
## selects variants with freq data from all VCF sources (1KG, ESP, gnomAD)
## also selects "older" rs numbers that are less likely to be retired and break these tests
# while(<>) { $line = $_; chomp($line); $num = $line =~ tr/\:/\:/; $rs = (split(qq{\t}, $line))[1]; next unless $rs =~ s/^rs//; print if $num > 14 && $rs < 1e7; }

## rand_lines.pl
# $num = shift; while(<>) { print if rand() < $num;}

## recut from http://www1.cuni.cz/~obo/textutils/
## a bit like cut but supports re-ordering and duplication of columns
##
  if($assembly eq 'GRCh38') {
    print OUT qq{1 1522770 1522770 T/C + rs1135025 C:0.0507 C:0.0504 C:0.0109 C:0.0109 C:0.0102 C:0.008872 C:0.007791 C:0.01547 C:0.006986 C:0.05807 C:0.007647 C:0.01034 C:0.01265 C:0.007168 C:0.01825 C:0.009816
1 3727810 3727810 A/G + rs2146657 G:0.7723 G:0.6758 G:0.3353 G:0.828 G:0.6074 G:0.7815 G:0.8001 G:0.6857 G:0.7607 G:0.6051 G:0.7937 G:0.3155 G:0.7581 G:0.7925 G:0.744 G:0.6427
1 11529854 11529854 G/A + rs2076468 A:0.0726 A:0.0994 A:0.4385 A:0.1064 A:0.2035 A:0.08973 A:0.1407 A:0.1544 A:0.09083 A:0.0864 A:0.1602 A:0.4474 A:0.1258 A:0.1351 A:0.1368 A:0.1877
1 12861684 12861684 T/G + rs1063798 G:0.2595 G:0.1787 G:0.1716 G:0.1799 G:0.1483 G:0.07423 G:0.0431 G:0.1391 G:0.2285 G:0.1752 G:0.1337 G:0.1998 G:0.09008 G:0.1187 G:0.1784 G:0.1291
1 16053748 16053748 A/G + rs2275166 G:0.82 G:0.5519 G:0.8056 G:0.6262 G:0.6186 G:0.7994 G:0.6503 G:0.6527 G:0.7984 G:0.5168 G:0.5828 G:0.7943 G:0.7409 G:0.6598 G:0.6518 G:0.5829
1 16208235 16208235 G/A + rs221059 A:0.1664 A:0.4063 A:0.1528 A:0.2654 A:0.4223 A:0.1754 A:0.261 A:0.2846 A:0.1716 A:0.455 A:0.3086 A:0.1264 A:0.1816 A:0.2581 A:0.2947 A:0.3993
1 21227230 21227230 G/T + rs3026903 T:0.0076 T:0.0173 T:0 T:0.0408 T:0.0102 T:0.01384 T:0.0314 T:0.02817 T:0.01065 T:0.01286 T:0.03594 T:0 T:0.06599 T:0.0355 T:0.02971 T:0.01257
1 21468895 21468895 A/G + rs1827293 G:0.1445 G:0.4755 G:0.5804 G:0.5417 G:0.4826 G:0.2202 G:0.553 G:0.5268 G:0.1962 G:0.5554 G:0.5769 G:0.5682 G:0.543 G:0.5579 G:0.5603 G:0.4902
1 26170396 26170396 C/T + rs2232655 T:0.0166 T:0.0014 T:0 T:0 T:0 T:0.009532 T:0.0001163 T:0.0007359 T:0.009022 T:0.0007459 T:0 T:0 T:0 T:0.0001082 T:0.0009181 T:0
1 27955695 27955695 C/T + rs3813804 T:0.1641 T:0.3112 T:0.253 T:0.2704 T:0.2485 T:0.2022 T:0.281 T:0.2709 T:0.1887 T:0.3273 T:0.2996 T:0.2448 T:0.216 T:0.2807 T:0.2592 T:0.2619
10 1048969 1048969 C/G + rs4880760 G:0.0015 G:0.0187 G:0.002 G:0.0209 G:0.0082 G:0.00379 G:0.02554 G:0.02285 G:0.004408 G:0.02059 G:0.08034 G:0 G:0.002882 G:0.03189 G:0.03601 G:0.006435
10 12816373 12816373 T/C + rs1644417 C:0.5514 C:0.1744 C:0.4077 C:0.0388 C:0.1268 C:0.4837 C:0.02977 C:0.1231 C:0.5105 C:0.1975 C:0.07517 C:0.4254 C:0.05358 C:0.02603 C:0.09952 C:0.1065
10 13204190 13204190 C/T + rs7072097 T:0.0711 T:0.0144 T:0.0437 T:0.0268 T:0.1033 T:0.07059 T:0.02767 T:0.03438 T:0.0773 T:0.01456 T:0.03317 T:0.04228 T:0.00677 T:0.02713 T:0.03011 T:0.07783
10 16694978 16694978 T/C + rs3740173 C:0.2474 C:0.415 C:0.4196 C:0.3996 C:0.2975 C:0.2837 C:0.3905 C:0.3738 C:0.2708 C:0.3807 C:0.4573 C:0.4385 C:0.3402 C:0.3908 C:0.3926 C:0.3179
10 16904126 16904126 T/C + rs7088988 C:0.0038 C:0.0288 C:0 C:0.0646 C:0.0511 C:0.01089 C:0.07779 C:0.05689 C:0.0117 C:0.02908 C:0.06792 C:0.0004058 C:0.04791 C:0.07717 C:0.05732 C:0.07067
10 26219214 26219214 C/A + rs2839672 A:0 A:0.0043 A:0 A:0.008 A:0.0031 A:0.001816 A:0.007442 A:0.005505 A:0.001242 A:0.001193 A:0.002742 A:0 A:0.0101 A:0.008214 A:0.004932 A:0.003226
10 35569041 35569041 T/C + rs7092436 C:0.0825 C:0.0058 C:0 C:0 C:0 C:0.07762 C:0.0002326 C:0.00564 C:0.07777 C:0.004498 C:0 C:0 C:0 C:0.0001639 C:0.002195 C:0.0001625
10 48242400 48242400 A/G + rs280612 G:0.6536 G:0.938 G:0.999 G:0.9453 G:0.9673 G:0.7253 G:0.9524 G:0.9365 G:0.6821 G:0.9611 G:0.9074 G:1 G:0.9543 G:0.9491 G:0.9352 G:0.9524
10 49493091 49493091 A/G + rs4253133 G:0.0673 G:0.0144 G:0 G:0 G:0 G:0.04925 G:0.0004651 G:0.003877 G:0.05024 G:0.0034 G:0.003155 G:0 G:0 G:0.0002252 G:0.002384 G:6.525e-05
10 60074398 60074398 G/A + rs3134609 A:0 A:0.0043 A:0.001 A:0.008 A:0.0153 A:0.001816 A:0.01465 A:0.009081 A:0.002418 A:0.003694 A:0.0001017 A:0.000116 A:0.003229 A:0.0139 A:0.00749 A:0.01326
11 244129 244129 C/T + rs7116130 T:0.8918 T:0.83 T:0.8036 T:0.7306 T:0.8006 T:0.8726 T:0.7284 T:0.7836 T:0.8913 T:0.8506 T:0.7099 T:0.8356 T:0.7709 T:0.7482 T:0.7518 T:0.7951
11 1836013 1836013 G/A + rs907608 A:0.0287 A:0.1758 A:0.0913 A:0.1392 A:0.2096 A:0.04678 A:0.1362 A:0.1532 A:0.04215 A:0.2094 A:0.1124 A:0.06896 A:0.2286 A:0.1427 A:0.1739 A:0.1931
11 2160036 2160036 C/T + rs3842750 T:0.0008 T:0 T:0 T:0 T:0 T:0.0009259 T:0 T:0.0001311 T:0.0003398 T:0 T:0 T:0 T:0 T:0 T:0 T:0.0007152
11 3038164 3038164 G/A + rs3205318 A:0.0303 A:0.0029 A:0 A:0 A:0 A:0.02816 A:0.0002327 A:0.002002 A:0.0264 A:0.001638 A:0.0008122 A:0 A:0 A:0.0001432 A:0.001641 A:3.249e-05
11 6627193 6627193 G/A + rs4758443 A:0.1233 A:0.3372 A:0.4236 A:0.3728 A:0.2975 A:0.1722 A:0.3678 A:0.3558 A:0.1484 A:0.345 A:0.4025 A:0.4437 A:0.396 A:0.377 A:0.3649 A:0.2989
11 8684859 8684859 T/G + rs4670 G:0.2617 G:0.1009 G:0.0635 G:0.0646 G:0.1339 G:0.2099 G:0.08566 G:0.1044 G:0.2245 G:0.11 G:0.1658 G:0.08894 G:0.0692 G:0.08281 G:0.1041 G:0.1312
11 18173280 18173280 C/G + rs2468774 G:0.2549 G:0.3415 G:0.4395 G:0.3062 G:0.3364 G:0.2456 G:0.2852 G:0.3059 G:0.2437 G:0.3433 G:0.2816 G:0.432 G:0.2414 G:0.2875 G:0.2839 G:0.35
11 22343198 22343198 G/A + rs2246710 A:0.2352 A:0.1571 A:0.0149 A:0.2406 A:0.1575 A:0.2267 A:0.22 A:0.1843 A:0.2301 A:0.1159 A:0.1474 A:0.01641 A:0.2091 A:0.2269 A:0.2049 A:0.1656
11 33559760 33559760 G/A + rs2281380 A:0.469 A:0.3213 A:0.3264 A:0.3887 A:0.3937 A:0.4592 A:0.386 A:0.3873 A:0.4718 A:0.3239 A:0.4432 A:0.3308 A:0.443 A:0.3807 A:0.4004 A:0.4097
11 57236107 57236107 C/T + rs7943508 T:0.0492 T:0.0086 T:0 T:0 T:0 T:0.04089 T:0.0003492 T:0.003224 T:0.04751 T:0.001579 T:0 T:0 T:0 T:2.691e-05 T:0.002007 T:3.249e-05
12 26715389 26715389 T/C + rs1900941 C:0.0295 C:0.2435 C:0.1944 C:0.492 C:0.1912 C:0.08249 C:0.4702 C:0.3523 C:0.07336 C:0.2209 C:0.3331 C:0.1834 C:0.4282 C:0.4758 C:0.3907 C:0.2257
12 27667369 27667369 C/T + rs2033024 T:0.938 T:0.8876 T:0.881 T:0.9533 T:0.9785 T:0.9274 T:0.9431 T:0.934 T:0.9339 T:0.8695 T:0.9428 T:0.8947 T:0.9449 T:0.9433 T:0.9302 T:0.9702
12 31132235 31132235 A/T + rs2287449 T:0.5356 T:0.6585 T:0.5655 T:0.7177 T:0.7055 T:0.5344 T:0.7261 T:0.6889 T:0.5187 T:0.6729 T:0.7113 T:0.5621 T:0.7437 T:0.7294 T:0.6971 T:0.6911
12 45016292 45016292 A/G + rs1224442 G:0.7065 G:0.732 G:0.7887 G:0.5825 G:0.8732 G:0.6675 G:0.5926 G:0.6699 G:0.6796 G:0.8117 G:0.5743 G:0.8192 G:0.6092 G:0.5905 G:0.6602 G:0.843
12 45350619 45350619 C/T + rs4108249 T:0.1293 T:0.0173 T:0 T:0 T:0.0092 T:0.1112 T:0.001628 T:0.0094 T:0.1131 T:0.00564 T:0.008319 T:5.969e-05 T:0 T:0.0008549 T:0.006401 T:0.007186
12 47077656 47077656 G/A + rs2269828 A:0.1346 A:0.2882 A:0.124 A:0.3211 A:0.2904 A:0.1716 A:0.3242 A:0.2771 A:0.1581 A:0.225 A:0.2754 A:0.08163 A:0.3134 A:0.3238 A:0.3062 A:0.3022
12 50350311 50350311 G/A + rs7296277 A:0.8888 A:0.9928 A:1 A:0.999 A:1 A:0.8909 A:1 A:0.9926 A:0.8972 A:0.9941 A:1 A:1 A:1 A:0.9998 A:0.994 A:0.9999
12 51242475 51242475 C/T + rs1049467 T:0.7905 T:0.8617 T:0.8988 T:0.7614 T:0.7751 T:0.798 T:0.7683 T:0.7819 T:0.7886 T:0.8547 T:0.773 T:0.8882 T:0.7363 T:0.7592 T:0.7874 T:0.7565
13 27963400 27963400 T/G + rs2504212 G:0.7542 G:0.7507 G:0.8145 G:0.8171 G:0.8998 G:0.7526 G:0.8195 G:0.823 G:0.7605 G:0.8008 G:0.8265 G:0.8051 G:0.85 G:0.8169 G:0.8163 G:0.8915
13 44434629 44434629 T/C + rs2234260 C:0.0008 C:0.0086 C:0.001 C:0.0149 C:0.0112 C:0.002508 C:0.0149 C:0.01038 C:0.002072 C:0.004305 C:0.003085 C:0 C:0.01941 C:0.01372 C:0.007754 C:0.00921
13 45338825 45338825 C/T + rs2234225 T:0.0023 T:0 T:0 T:0 T:0.001 T:0.002951 T:0 T:0.000534 T:0.003315 T:0.0001158 T:0 T:0 T:0 T:0 T:0 T:0.00286
13 45340193 45340193 A/G + rs2234220 G:0.0242 G:0.0014 G:0 G:0 G:0 G:0.0168 G:0 G:0.001167 G:0.01606 G:0.0007166 G:0 G:0 G:0 G:3.641e-05 G:0.001687 G:3.363e-05
13 52712815 52712815 C/T + rs7330220 T:0.1452 T:0.0562 T:0 T:0.0726 T:0.047 T:0.1409 T:0.08477 T:0.06832 T:0.1352 T:0.03328 T:0.08995 T:0.0003639 T:0.04835 T:0.08398 T:0.07081 T:0.05668
13 98474798 98474798 C/T + rs9582232 T:0.0961 T:0.1427 T:0.126 T:0.1839 T:0.1493 T:0.09986 T:0.1887 T:0.167 T:0.09516 T:0.1198 T:0.2264 T:0.139 T:0.1683 T:0.1931 T:0.1743 T:0.1578
13 101454113 101454113 C/T + rs2297701 T:0.2201 T:0.2305 T:0.3839 T:0.1312 T:0.2219 T:0.199 T:0.1365 T:0.1899 T:0.1959 T:0.3099 T:0.2184 T:0.3823 T:0.1206 T:0.1314 T:0.1848 T:0.1897
13 102798700 102798700 C/A + rs7993350 A:0.1513 A:0.0231 A:0.1121 A:0.0507 A:0.0491 A:0.1328 A:0.03837 A:0.05221 A:0.1348 A:0.02559 A:0.04193 A:0.1012 A:0.04999 A:0.04218 A:0.04556 A:0.05004
13 113088281 113088281 C/T + rs9324220 T:0.2511 T:0.1182 T:0.004 T:0.1531 T:0.2362 T:0.2138 T:0.1529 T:0.1425 T:0.2093 T:0.0838 T:0.2121 T:0.00429 T:0.07327 T:0.1528 T:0.1492 T:0.2401
13 113118777 113118777 C/T + rs6044 T:0.0068 T:0.0014 T:0 T:0 T:0 T:0.007036 T:0 T:0.0005021 T:0.006867 T:0.0001491 T:0 T:0 T:0 T:4.508e-05 T:0.0001836 T:0.0002603
14 23080576 23080576 T/C + rs3811182 C:0.6914 C:0.3689 C:0.4554 C:0.3777 C:0.5102 C:0.6512 C:0.412 C:0.4323 C:0.6601 C:0.3541 C:0.4261 C:0.4271 C:0.3503 C:0.4166 C:0.4112 C:0.529
14 24050987 24050987 C/T + rs4982850 T:0 T:0.0014 T:0 T:0.002 T:0 T:0.0004539 T:0.002907 T:0.001801 T:0.0001978 T:0.002476 T:0 T:0 T:0.0002252 T:0.002856 T:0.004392 T:0.0002941
14 24262068 24262068 G/A + rs2228336 A:0.0242 A:0.0029 A:0 A:0 A:0 A:0.02134 A:0.0001163 A:0.001384 A:0.02048 A:0.0005361 A:0 A:0 A:0 A:2.727e-05 A:0.000548 A:9.747e-05
14 24437328 24437328 G/A + rs3742522 A:0.3933 A:0.1138 A:0.125 A:0.0527 A:0.0665 A:0.333 A:0.04085 A:0.07964 A:0.3584 A:0.07653 A:0.0449 A:0.13 A:0.03847 A:0.04439 A:0.06453 A:0.0705
14 51249470 51249470 T/G + rs7161242 G:0.6536 G:0.7262 G:0.6349 G:0.7495 G:0.6564 G:0.6595 G:0.7709 G:0.7351 G:0.6468 G:0.7374 G:0.8049 G:0.6681 G:0.7484 G:0.771 G:0.737 G:0.64
14 56605620 56605620 C/T + rs3737170 T:0.2466 T:0.0591 T:0.1825 T:0.0398 T:0.0532 T:0.2312 T:0.04852 T:0.06816 T:0.2245 T:0.0606 T:0.09911 T:0.1354 T:0.04942 T:0.04477 T:0.06067 T:0.05417
14 60144930 60144930 C/T + rs399535 T:0.5893 T:0.1571 T:0.3472 T:0.1511 T:0.4346 T:0.5161 T:0.1527 T:0.2107 T:0.5204 T:0.1187 T:0.2236 T:0.3012 T:0.1096 T:0.159 T:0.1916 T:0.3777
14 64468652 64468652 A/C + rs2230490 C:0.1929 C:0.0115 C:0.006 C:0.0089 C:0 C:0.1521 C:0.0003488 C:0.0111 C:0.159 C:0.006884 C:0.0006101 C:0 C:0 C:0.0002872 C:0.005846 C:0.0002601
14 65616075 65616075 A/C + rs2229678 C:0 C:0.0029 C:0 C:0.0268 C:0.001 C:0.00227 C:0.01349 C:0.01485 C:0.002028 C:0.006804 C:0.009957 C:0 C:0.07195 C:0.01407 C:0.01389 C:0.001496
14 71588751 71588751 A/C + rs1859643 C:0.9622 C:0.8473 C:0.9355 C:0.7614 C:0.8742 C:0.924 C:0.7547 C:0.8172 C:0.936 C:0.8782 C:0.7234 C:0.9401 C:0.8325 C:0.7564 C:0.7834 C:0.8671
15 40564790 40564790 C/T + rs3803354 T:0.7458 T:0.9121 T:0.872 T:0.9304 T:0.8722 T:0.7809 T:0.9301 T:0.9109 T:0.768 T:0.9253 T:0.9623 T:0.8741 T:0.9256 T:0.9337 T:0.9343 T:0.8731
15 41736156 41736156 C/A + rs2577960 A:0.8805 A:0.5519 A:0.8125 A:0.6193 A:0.6544 A:0.8301 A:0.6458 A:0.6501 A:0.8265 A:0.5813 A:0.6316 A:0.796 A:0.5836 A:0.626 A:0.6199 A:0.6847
15 41847444 41847444 C/T + rs883329 T:0.0477 T:0.1873 T:0.6587 T:0.1243 T:0.1227 T:0.06605 T:0.142 T:0.1887 T:0.0673 T:0.2707 T:0.1444 T:0.6568 T:0.1304 T:0.1388 T:0.1675 T:0.1342
15 43017231 43017231 G/A + rs7182063 A:0.0348 A:0.0058 A:0 A:0.002 A:0.001 A:0.03638 A:0.002094 A:0.003765 A:0.03762 A:0.00266 A:0.004583 A:0 A:0.0001224 A:0.001212 A:0.004488 A:0.001487
15 52611575 52611575 T/C + rs8036680 C:0.1377 C:0.0706 C:0.001 C:0.0865 C:0.0736 C:0.1259 C:0.0981 C:0.08873 C:0.1321 C:0.06739 C:0.2172 C:0.0004058 C:0.08991 C:0.09352 C:0.1103 C:0.07686
15 62164746 62164746 G/A + rs8025811 A:0.3419 A:0.1196 A:0.001 A:0.1481 A:0.1462 A:0.2457 A:0.1144 A:0.1367 A:0.3279 A:0.08498 A:0.1642 A:0.0008772 A:0.1121 A:0.1518 A:0.1639 A:0.1624
15 63971585 63971585 G/A + rs7180050 A:0.18 A:0.415 A:0.0109 A:0.4374 A:0.2362 A:0.2347 A:0.4331 A:0.3458 A:0.2162 A:0.3362 A:0.4013 A:0.01109 A:0.4041 A:0.4287 A:0.3926 A:0.245
15 76883859 76883859 A/T + rs3812909 T:0.3427 T:0.3228 T:0.5774 T:0.3101 T:0.5286 T:0.3158 T:0.2849 T:0.3025 T:0.2705 T:0.2944 T:0.3406 T:0.4414 T:0.2987 T:0.2551 T:0.2766 T:0.4484
15 78925289 78925289 A/G + rs2289695 G:0.1036 G:0.1081 G:0.2817 G:0.1004 G:0.0665 G:0.09768 G:0.09655 G:0.1075 G:0.09902 G:0.09477 G:0.0947 G:0.2865 G:0.0826 G:0.09822 G:0.1048 G:0.08116
15 88878096 88878096 A/G + rs2280463 G:0.8139 G:0.3228 G:0.256 G:0.2932 G:0.3476 G:0.7475 G:0.31 G:0.3305 G:0.7636 G:0.3171 G:0.3154 G:0.2363 G:0.241 G:0.301 G:0.3329 G:0.3583
16 62777 62777 G/A + rs3213511 A:0 A:0 A:0.0159 A:0 A:0 A:0 A:0.0001163 A:0.001044 A:0 A:0 A:0 A:0.0145 A:0 A:8.954e-06 A:0.000547 A:9.747e-05
16 505505 505505 G/T + rs7201440 T:0.1135 T:0.0346 T:0.0079 T:0.0318 T:0.0051 T:0.1091 T:0.02846 T:0.02908 T:0.104 T:0.0257 T:0.02245 T:0.004604 T:0.03797 T:0.02863 T:0.02728 T:0.009285
16 879711 879711 C/T + rs2076425 T:0.208 T:0.0735 T:0.3026 T:0.1113 T:0.1472 T:0.1843 T:0.1036 T:0.1177 T:0.2096 T:0.05375 T:0.07417 T:0.2594 T:0.07411 T:0.1151 T:0.09943 T:0.1263
16 2770881 2770881 G/T + rs1050134 T:0.2632 T:0.1326 T:0.2202 T:0.166 T:0.092 T:0.2455 T:0.1519 T:0.157 T:0.247 T:0.1205 T:0.2587 T:0.2494 T:0.1587 T:0.1511 T:0.1597 T:0.08746
16 3021820 3021820 C/A + rs2232799 A:0.0545 A:0.0043 A:0.001 A:0.002 A:0 A:0.03913 A:0.0001163 A:0.002879 A:0.04022 A:0.001735 A:0 A:5.805e-05 A:0 A:0.0002068 A:0.001653 A:3.257e-05
16 3440922 3440922 G/C + rs2270494 C:0.1354 C:0.6066 C:0.6369 C:0.5547 C:0.5419 C:0.2014 C:0.5453 C:0.5446 C:0.1814 C:0.6683 C:0.4595 C:0.6416 C:0.5932 C:0.5507 C:0.5462 C:0.5053
16 8635096 8635096 C/T + rs1731044 T:0.9493 T:0.9928 T:1 T:1 T:1 T:0.9596 T:1 T:0.9971 T:0.957 T:0.9987 T:1 T:1 T:1 T:1 T:0.9987 T:1
16 11676107 11676107 C/T + rs1050069 T:0.0688 T:0.2032 T:0.252 T:0.4672 T:0.3323 T:0.129 T:0.469 T:0.3842 T:0.1102 T:0.1903 T:0.5226 T:0.2401 T:0.4676 T:0.4828 T:0.4065 T:0.3476
16 20429731 20429731 C/G + rs8062344 G:0.6846 G:0.2046 G:0.5774 G:0.2813 G:0.2935 G:0.6099 G:0.2945 G:0.3145 G:0.6361 G:0.1888 G:0.2287 G:0.5866 G:0.2995 G:0.2949 G:0.2818 G:0.2546
16 20784967 20784967 T/C + rs8058407 C:0.0166 C:0 C:0 C:0 C:0 C:0.01567 C:0 C:0.0008509 C:0.01271 C:0.0002906 C:0 C:0 C:0 C:0 C:0.000189 C:3.404e-05
17 410508 410508 G/A + rs4581766 A:0.2905 A:0.1657 A:0.2073 A:0.1093 A:0.1779 A:0.1854 A:0.08375 A:0.1472 A:0.2667 A:0.1629 A:0.2212 A:0.2485 A:0.1112 A:0.1014 A:0.1649 A:0.1362
17 1006211 1006211 A/T + rs2586306 T:0.4198 T:0.4957 T:0.2946 T:0.506 T:0.3691 T:0.4653 T:0.5391 T:0.4721 T:0.4441 T:0.4539 T:0.5702 T:0.3016 T:0.497 T:0.5197 T:0.5086 T:0.4007
17 2415401 2415401 G/A + rs4790335 A:0.3752 A:0.7666 A:0.9643 A:0.6571 A:0.8282 A:0.4417 A:0.6721 A:0.7205 A:0.4111 A:0.8268 A:0.702 A:0.9767 A:0.7361 A:0.6646 A:0.7016 A:0.826
17 3724325 3724325 C/T + rs1185511 T:0.587 T:0.0403 T:0 T:0 T:0 T:0.4255 T:0.002998 T:0.03421 T:0.497 T:0.02467 T:0.00588 T:0.0001373 T:8.218e-05 T:0.001675 T:0.02186 T:0.001529
17 3925408 3925408 T/C + rs887387 C:0.1823 C:0.4424 C:0.2857 C:0.335 C:0.2669 C:0.222 C:0.3799 C:0.3526 C:0.2048 C:0.4449 C:0.3945 C:0.2831 C:0.322 C:0.3766 C:0.3475 C:0.2877
17 4077854 4077854 G/C + rs8074454 C:0.4289 C:0.2536 C:0.1835 C:0.3221 C:0.2239 C:0.3976 C:0.3291 C:0.2925 C:0.4171 C:0.2508 C:0.2637 C:0.1717 C:0.2994 C:0.3194 C:0.3213 C:0.2461
17 4539780 4539780 T/A + rs9905742 A:0.003 A:0.0317 A:0 A:0.0537 A:0.0266 A:0.00522 A:0.04023 A:0.03043 A:0.004641 A:0.02567 A:0.08012 A:5.798e-05 A:0.01253 A:0.03989 A:0.0398 A:0.02602
17 6635035 6635035 G/A + rs2301873 A:0.2655 A:0.3271 A:0.2897 A:0.4543 A:0.365 A:0.2902 A:0.427 A:0.3844 A:0.3019 A:0.258 A:0.4897 A:0.2994 A:0.4471 A:0.4304 A:0.3885 A:0.3637
17 7934938 7934938 T/C + rs7217064 C:0.0144 C:0 C:0 C:0 C:0 C:0.01754 C:0 C:0.001683 C:0.0185 C:0.001121 C:0 C:0 C:0 C:4.735e-05 C:0.0005656 C:0
17 9646124 9646124 C/T + rs3744756 T:0.1664 T:0.1124 T:0.4187 T:0.0149 T:0.3088 T:0.0825 T:0.002974 T:0.105 T:0.1727 T:0.1666 T:0.01948 T:0.4415 T:0.04838 T:0.004307 T:0.06977 T:0.258
18 3214920 3214920 G/C + rs7232679 C:0.388 C:0.2478 C:0.2252 C:0.2416 C:0.2843 C:0.3509 C:0.2615 C:0.2724 C:0.3776 C:0.2124 C:0.2158 C:0.2057 C:0.3056 C:0.2774 C:0.2635 C:0.2993
18 11752944 11752944 C/A + rs1895689 A:0.5318 A:0.1614 A:0.3284 A:0.2664 A:0.3364 A:0.4753 A:0.2878 A:0.281 A:0.4927 A:0.1428 A:0.3377 A:0.3484 A:0.1746 A:0.2811 A:0.2779 A:0.349
18 32310334 32310334 G/C + rs683701 C:0.3026 C:0.1326 C:0.0486 C:0.1262 C:0.1902 C:0.3046 C:0.1474 C:0.143 C:0.3106 C:0.08813 C:0.1327 C:0.03819 C:0.1471 C:0.1451 C:0.1293 C:0.1737
18 33744970 33744970 A/G + rs7232237 G:0.7352 G:0.6816 G:0.8988 G:0.5189 G:0.5399 G:0.7378 G:0.5134 G:0.5925 G:0.7356 G:0.7592 G:0.5535 G:0.9103 G:0.5335 G:0.5022 G:0.5703 G:0.5481
18 35977503 35977503 A/G + rs2276314 G:0.323 G:0.1859 G:0.248 G:0.2137 G:0.2352 G:0.3155 G:0.2229 G:0.2286 G:0.3197 G:0.1736 G:0.2857 G:0.253 G:0.2216 G:0.2231 G:0.2159 G:0.2382
18 45946806 45946806 T/C + rs8090934 C:0.1309 C:0.0231 C:0 C:0.003 C:0.0031 C:0.1055 C:0.0018 C:0.009794 C:0.1128 C:0.01012 C:0.002558 C:0 C:0 C:0.001979 C:0.008984 C:0.00153
18 54293619 54293619 G/C + rs660793 C:0.2186 C:0.1988 C:0.0347 C:0.2734 C:0.2587 C:0.2286 C:0.2758 C:0.2368 C:0.2389 C:0.1691 C:0.3122 C:0.02198 C:0.2331 C:0.288 C:0.262 C:0.2526
18 58481816 58481816 A/G + rs3744867 G:0.4372 G:0.3343 G:0.2143 G:0.4692 G:0.3231 G:0.414 G:0.423 G:0.3953 G:0.4248 G:0.2931 G:0.3916 G:0.2141 G:0.5252 G:0.4363 G:0.4199 G:0.3478
18 58537759 58537759 C/T + rs3809970 T:0.4312 T:0.4265 T:0.2361 T:0.5716 T:0.3497 T:0.4501 T:0.5337 T:0.4656 T:0.451 T:0.3629 T:0.425 T:0.2229 T:0.5392 T:0.5498 T:0.4927 T:0.3709
18 76904418 76904418 A/G + rs3794873 G:0.0008 G:0.0144 G:0.0288 G:0.0229 G:0.0174 G:0.006343 G:0.02834 G:0.02406 G:0.004152 G:0.01121 G:0.03439 G:0.02797 G:0.02827 G:0.02722 G:0.02435 G:0.02736
19 1041853 1041853 G/T + rs3764644 T:0.1513 T:0.0187 T:0.0923 T:0.0358 T:0.0634 T:0.1233 T:0.03907 T:0.04789 T:0.1356 T:0.01946 T:0.02782 T:0.08507 T:0.02679 T:0.03998 T:0.03315 T:0.06529
19 1555552 1555552 T/C + rs2277752 C:0.0537 C:0.1311 C:0.3194 C:0.1173 C:0.1258 C:0.06529 C:0.09532 C:0.1306 C:0.07058 C:0.114 C:0.1351 C:0.3318 C:0.161 C:0.1081 C:0.1066 C:0.1169
19 3586473 3586473 A/G + rs8100350 G:0.6293 G:0.3285 G:0.2073 G:0.2545 G:0.3415 G:0.5824 G:0.2411 G:0.2851 G:0.5975 G:0.3236 G:0.1994 G:0.2316 G:0.277 G:0.2376 G:0.26 G:0.3286
19 4054210 4054210 C/T + rs7251080 T:0.1876 T:0.0187 T:0 T:0 T:0 T:0.1791 T:0.0005814 T:0.0135 T:0.1864 T:0.007781 T:0 T:0 T:0 T:0.0003379 T:0.006312 T:0.0003744
19 4288223 4288223 C/T + rs888929 T:0.618 T:0.5591 T:0.3393 T:0.3777 T:0.5041 T:0.5733 T:0.3667 T:0.4131 T:0.5794 T:0.5685 T:0.363 T:0.3346 T:0.3177 T:0.3631 T:0.3897 T:0.4767
19 4493996 4493996 C/T + rs2288931 T:0 T:0 T:0.0069 T:0 T:0 T:0 T:0.0001215 T:0.0001223 T:0 T:0.0003016 T:0 T:0.00106 T:0 T:9.407e-06 T:0 T:0
19 5222820 5222820 G/A + rs2230610 A:0 A:0.0058 A:0 A:0.0298 A:0.001 A:0.004344 A:0.02789 A:0.01938 A:0.004779 A:0.007583 A:0.008478 A:0 A:0.03097 A:0.03349 A:0.02098 A:0.001453
19 5610034 5610034 G/A + rs806706 A:0.0015 A:0.0231 A:0 A:0.0298 A:0.0123 A:0.006582 A:0.02872 A:0.02531 A:0.005566 A:0.0137 A:0.03557 A:0 A:0.04886 A:0.03231 A:0.02627 A:0.01608
19 6502126 6502126 C/A + rs390530 A:0.0968 A:0.0519 A:0.2381 A:0.0586 A:0.1861 A:0.08203 A:0.05081 A:0.09715 A:0.09375 A:0.06754 A:0.08404 A:0.2235 A:0.102 A:0.05459 A:0.07781 A:0.2057
19 7847930 7847930 C/T + rs2307003 T:0.1672 T:0.5159 T:0.4583 T:0.3618 T:0.7157 T:0.2024 T:0.3457 T:0.4507 T:0.2072 T:0.5339 T:0.3932 T:0.4417 T:0.4887 T:0.3681 T:0.4211 T:0.6528
2 276942 276942 A/G + rs7573495 G:0.1823 G:0.2565 G:0.5853 G:0.3559 G:0.2454 G:0.1996 G:0.339 G:0.3338 G:0.2011 G:0.2439 G:0.3351 G:0.5964 G:0.3762 G:0.3428 G:0.3335 G:0.2721
2 3613266 3613266 G/A + rs7567724 A:0.1051 A:0.0663 A:0.0228 A:0.1103 A:0.089 A:0.09973 A:0.1205 A:0.1098 A:0.1038 A:0.05878 A:0.2287 A:0.02494 A:0.1234 A:0.1244 A:0.1185 A:0.1153
2 10607145 10607145 A/G + rs6759740 G:0.7738 G:0.5634 G:0.624 G:0.5169 G:0.3988 G:0.7578 G:0.5148 G:0.5391 G:0.7554 G:0.5681 G:0.3736 G:0.6273 G:0.6512 G:0.5142 G:0.499 G:0.4257
2 10913614 10913614 C/T + rs3732105 T:0.3389 T:0.4928 T:0.5099 T:0.6322 T:0.4898 T:0.3761 T:0.6398 T:0.5764 T:0.3588 T:0.4567 T:0.6508 T:0.5331 T:0.5836 T:0.6507 T:0.5856 T:0.5398
2 17781251 17781251 C/T + rs300169 T:0.5719 T:0.7608 T:0.9851 T:0.6372 T:0.8978 T:0.6074 T:0.6245 T:0.7038 T:0.5955 T:0.828 T:0.6545 T:0.9844 T:0.5719 T:0.6225 T:0.6818 T:0.8732
2 26582893 26582893 C/T + rs7595964 T:0.1483 T:0.0231 T:0 T:0.001 T:0.0245 T:0.128 T:0.002791 T:0.0151 T:0.1366 T:0.008669 T:0.008731 T:5.798e-05 T:8.969e-05 T:0.002158 T:0.01222 T:0.03056
2 28538415 28538415 A/G + rs7589254 G:0.6619 G:0.7853 G:0.9931 G:0.9254 G:0.9417 G:0.721 G:0.9275 G:0.8976 G:0.6873 G:0.7635 G:0.8417 G:0.9942 G:0.979 G:0.9343 G:0.9081 G:0.9234
2 31177173 31177173 A/G + rs4516476 G:0.6225 G:0.4481 G:0.1756 G:0.2256 G:0.2894 G:0.5434 G:0.2231 G:0.2914 G:0.5729 G:0.4488 G:0.2071 G:0.1758 G:0.2646 G:0.2269 G:0.2677 G:0.2958
2 31201827 31201827 C/A + rs749251 A:0.059 A:0.0101 A:0.0228 A:0.0338 A:0.0041 A:0.05202 A:0.03426 A:0.02649 A:0.05823 A:0.01454 A:0.0021 A:0.01312 A:0.03818 A:0.03586 A:0.02405 A:0.01194
2 31403179 31403179 G/C + rs1346644 C:0.1293 C:0.1081 C:0.0625 C:0.168 C:0.18 C:0.1414 C:0.1424 C:0.1272 C:0.147 C:0.07629 C:0.1247 C:0.04832 C:0.1233 C:0.1448 C:0.1299 C:0.1564
20 2580070 2580070 T/G + rs2422794 G:0.1536 G:0.036 G:0.0288 G:0.0994 G:0.047 G:0.1532 G:0.101 G:0.08409 G:0.155 G:0.03466 G:0.07712 G:0.02074 G:0.117 G:0.1017 G:0.08496 G:0.05157
20 3147757 3147757 T/C + rs3746698 C:0.0295 C:0.0029 C:0.0099 C:0.008 C:0.0634 C:0.02678 C:0.005233 C:0.01322 C:0.02483 C:0.003842 C:0.01604 C:0.01194 C:0.0009417 C:0.005456 C:0.00802 C:0.05549
20 9529616 9529616 T/A + rs2232268 A:0.7035 A:0.2089 A:0.2123 A:0.1879 A:0.2076 A:0.6158 A:0.1928 A:0.2117 A:0.6351 A:0.1222 A:0.2371 A:0.2272 A:0.1707 A:0.1847 A:0.1935 A:0.2104
20 9563041 9563041 T/G + rs2297347 G:0.8389 G:0.2219 G:0.245 G:0.2286 G:0.2362 G:0.7392 G:0.2248 G:0.2476 G:0.7646 G:0.1369 G:0.2963 G:0.2771 G:0.2014 G:0.2156 G:0.2295 G:0.2317
20 10038445 10038445 A/G + rs575534 G:0.7526 G:0.5245 G:0.3145 G:0.6024 G:0.5624 G:0.7224 G:0.6051 G:0.5697 G:0.7346 G:0.5433 G:0.515 G:0.3131 G:0.5888 G:0.5995 G:0.5766 G:0.5547
20 20228327 20228327 G/A + rs6046740 A:0.09 A:0.0173 A:0 A:0 A:0 A:0.06537 A:0.0005814 A:0.005354 A:0.06797 A:0.00539 A:0.0002031 A:0 A:0 A:0.0005015 A:0.004927 A:0.0003899
20 23036333 23036333 T/G + rs3746726 G:0.2179 G:0.281 G:0.2351 G:0.3698 G:0.4479 G:0.234 G:0.3862 G:0.3517 G:0.2316 G:0.2865 G:0.4958 G:0.2419 G:0.3275 G:0.3719 G:0.3676 G:0.439
20 33169272 33169272 G/A + rs6059139 A:0.0998 A:0.0058 A:0 A:0.005 A:0 A:0.0926 A:0.001977 A:0.00835 A:0.1025 A:0.005986 A:0.003452 A:0 A:0 A:0.001764 A:0.007291 A:0.0004874
20 33710517 33710517 C/T + rs2626529 T:0.6838 T:0.5447 T:0.7887 T:0.4811 T:0.59 T:0.6332 T:0.4799 T:0.5315 T:0.6333 T:0.5923 T:0.4776 T:0.7993 T:0.4997 T:0.4594 T:0.5115 T:0.5481
20 34090562 34090562 T/G + rs4911406 G:0.0015 G:0.0014 G:0.0516 G:0.004 G:0.0031 G:0.0009217 G:0.006098 G:0.008044 G:0.0008768 G:0.002051 G:0.0009292 G:0.04143 G:0.01547 G:0.006111 G:0.005871 G:0.003967
21 26837711 26837711 A/G + rs3183693 G:0 G:0.013 G:0 G:0.0189 G:0 G:0.00227 G:0.01942 G:0.01122 G:0.003136 G:0.007594 G:0.007411 G:0 G:0.004439 G:0.01965 G:0.01003 G:0.001267
21 41275363 41275363 C/T + rs6517660 T:0.0575 T:0.0086 T:0 T:0 T:0 T:0.05107 T:0.0003488 T:0.004333 T:0.05528 T:0.00405 T:0 T:5.798e-05 T:0 T:0.0004126 T:0.002007 T:0.0008446
21 42119183 42119183 G/A + rs220146 A:0.1543 A:0.33 A:0.2222 A:0.2227 A:0.2209 A:0.1519 A:0.2116 A:0.2271 A:0.1539 A:0.3175 A:0.2583 A:0.22 A:0.2103 A:0.2171 A:0.2351 A:0.2061
21 42137588 42137588 A/C + rs3819142 C:0.1452 C:0.245 C:0.256 C:0.1769 C:0.3139 C:0.1566 C:0.1798 C:0.2083 C:0.1671 C:0.2219 C:0.1823 C:0.2386 C:0.2302 C:0.1806 C:0.1872 C:0.2939
21 44600636 44600636 T/C + rs512211 C:0.0741 C:0.0058 C:0.001 C:0 C:0 C:0.0535 C:0.0003538 C:0.004196 C:0.05932 C:0.002421 C:0 C:0.0006984 C:0 C:0.000117 C:0.002206 C:0.000196
21 44638033 44638033 G/A + rs7280841 A:0.0779 A:0.0086 A:0 A:0 A:0.0031 A:0.06764 A:0.0001163 A:0.00541 A:0.07162 A:0.003157 A:0 A:0.000116 A:0.0002242 A:0.0002063 A:0.003466 A:0.002599
22 17188309 17188309 G/A + rs7284498 A:0.1755 A:0.1282 A:0.0962 A:0.2346 A:0.2597 A:0.1614 A:0.2149 A:0.1953 A:0.1618 A:0.1016 A:0.3085 A:0.104 A:0.2067 A:0.2206 A:0.2119 A:0.2275
22 19137353 19137353 G/A + rs3765610 A:0.0015 A:0.0634 A:0.0685 A:0.0149 A:0.0757 A:0.003177 A:0.001744 A:0.03264 A:0.002027 A:0.1026 A:0.003156 A:0.05187 A:0.08148 A:0.002941 A:0.03069 A:0.04301
22 20108841 20108841 C/T + rs2286927 T:0.0008 T:0 T:0.004 T:0 T:0 T:0.0004539 T:0.0002326 T:0.0004068 T:0.0003935 T:8.943e-05 T:0.003571 T:0.002496 T:0 T:5.429e-05 T:0.0003668 T:0.00013
22 20974834 20974834 A/G + rs733455 G:0.357 G:0.6916 G:0.494 G:0.8469 G:0.6953 G:0.4262 G:0.8363 G:0.7511 G:0.4101 G:0.7019 G:0.8452 G:0.5108 G:0.8197 G:0.8442 G:0.8002 G:0.6868
22 21701105 21701105 A/G + rs1860 G:0.6861 G:0.6297 G:0.7411 G:0.6183 G:0.6544 G:0.6682 G:0.6034 G:0.6425 G:0.667 G:0.6391 G:0.6867 G:0.7516 G:0.7349 G:0.6046 G:0.6386 G:0.63
22 31714952 31714952 C/T + rs9619227 T:0.0598 T:0.0029 T:0 T:0 T:0 T:0.04697 T:0.0006285 T:0.003048 T:0.04373 T:0.002962 T:0 T:0 T:0 T:0.0005116 T:0.003801 T:0.0002191
22 32193036 32193036 C/T + rs136478 T:0.5227 T:0.2723 T:0.1419 T:0.4642 T:0.3599 T:0.498 T:0.4521 T:0.3871 T:0.5033 T:0.2429 T:0.4685 T:0.1386 T:0.3449 T:0.4613 T:0.407 T:0.3579
22 32493306 32493306 C/T + rs5998515 T:0.0106 T:0 T:0 T:0.002 T:0 T:0.007944 T:0.0004651 T:0.0007236 T:0.008894 T:0.0005064 T:0 T:0 T:0 T:0.0002062 T:0.0003651 T:0
22 32518246 32518246 C/T + rs5998527 T:0.3812 T:0.0331 T:0 T:0.008 T:0.0818 T:0.2798 T:0.002791 T:0.03186 T:0.3127 T:0.01902 T:0.006004 T:0.000116 T:0 T:0.00286 T:0.02417 T:0.0605
22 35412715 35412715 C/T + rs2307344 T:0.0832 T:0.0735 T:0.1458 T:0.0736 T:0.0429 T:0.08379 T:0.07103 T:0.07946 T:0.08485 T:0.07071 T:0.1161 T:0.1351 T:0.08608 T:0.07585 T:0.07849 T:0.04352
3 10092265 10092265 A/G + rs9811771 G:0.5431 G:0.1974 G:0.0675 G:0.1412 G:0.2127 G:0.4902 G:0.1653 G:0.1844 G:0.5032 G:0.1832 G:0.1593 G:0.08999 G:0.1211 G:0.1599 G:0.167 G:0.2259
3 10178045 10178045 A/G + rs1681663 G:0.6286 G:0.6916 G:0.2192 G:0.7694 G:0.6258 G:0.6548 G:0.7747 G:0.6966 G:0.6465 G:0.6788 G:0.7824 G:0.2536 G:0.7947 G:0.7724 G:0.7349 G:0.6335
3 10340464 10340464 G/A + rs154242 A:0.7443 A:0.7363 A:0.8849 A:0.5845 A:0.8425 A:0.6798 A:0.5815 A:0.6641 A:0.7054 A:0.7685 A:0.5746 A:0.8881 A:0.5752 A:0.5784 A:0.6241 A:0.8151
3 16434024 16434024 C/T + rs6442602 T:0.0061 T:0.0014 T:0 T:0 T:0 T:0.0134 T:0 T:0.001085 T:0.01518 T:0.000631 T:0 T:0 T:0 T:9.585e-06 T:0.000788 T:0
3 17089718 17089718 C/A + rs9822682 A:0.382 A:0.1037 A:0.0327 A:0.1004 A:0.1493 A:0.3094 A:0.1108 A:0.1177 A:0.3241 A:0.09032 A:0.08761 A:0.02883 A:0.1141 A:0.1074 A:0.1095 A:0.1452
3 41249590 41249590 G/A + rs2293304 A:0.0908 A:0.3055 A:0.2569 A:0.4085 A:0.3957 A:0.1424 A:0.4409 A:0.3791 A:0.1359 A:0.3386 A:0.3629 A:0.2483 A:0.4146 A:0.4393 A:0.3891 A:0.3793
3 43718570 43718570 T/G + rs887472 G:0.1475 G:0.0735 G:0.0982 G:0.0288 G:0.1472 G:0.1323 G:0.0243 G:0.05425 G:0.1359 G:0.07812 G:0.01931 G:0.1032 G:0.01359 G:0.02383 G:0.04881 G:0.1119
3 47928244 47928244 G/A + rs20565 A:0.329 A:0.0375 A:0 A:0.001 A:0 A:0.3057 A:0.002093 A:0.02259 A:0.3074 A:0.01721 A:0.004365 A:0 A:0 A:0.001271 A:0.01605 A:0.0001949
3 48660221 48660221 C/G + rs3821875 G:0.1271 G:0.0793 G:0.1171 G:0.1193 G:0.0409 G:0.1337 G:0.1179 G:0.09971 G:0.1269 G:0.05352 G:0.1264 G:0.08697 G:0.1321 G:0.1142 G:0.1032 G:0.05848
3 49687375 49687375 T/C + rs3020779 C:0.7095 C:0.8573 C:0.9871 C:0.834 C:0.9622 C:0.732 C:0.8161 C:0.8639 C:0.7279 C:0.9115 C:0.8509 C:0.9931 C:0.9087 C:0.819 C:0.8426 C:0.945
4 947730 947730 C/T + rs2290402 T:0.025 T:0.1167 T:0.3056 T:0.1521 T:0.182 T:0.0469 T:0.162 T:0.154 T:0.04159 T:0.1073 T:0.2651 T:0.3098 T:0.05217 T:0.1561 T:0.1625 T:0.2061
4 1728572 1728572 A/G + rs798759 G:0.3245 G:0.3184 G:0.1597 G:0.4394 G:0.2873 G:0.3545 G:0.4507 G:0.3721 G:0.3364 G:0.307 G:0.3935 G:0.18 G:0.3484 G:0.4469 G:0.3797 G:0.3061
4 7726745 7726745 C/T + rs2301746 T:0.497 T:0.3934 T:0.5893 T:0.5388 T:0.5552 T:0.4775 T:0.4918 T:0.4841 T:0.4877 T:0.3482 T:0.4982 T:0.5832 T:0.4851 T:0.4981 T:0.4845 T:0.519
4 7740173 7740173 G/A + rs3733602 A:0 A:0 A:0.0665 A:0 A:0.0031 A:0 A:0.0001192 A:0.006632 A:7.362e-05 A:0.0002497 A:0.0008414 A:0.08173 A:0.001211 A:0.000457 A:0.003492 A:0.002234
4 15715698 15715698 G/A + rs3213710 A:0.2315 A:0.5749 A:0.4266 A:0.5308 A:0.5644 A:0.2705 A:0.523 A:0.5187 A:0.2494 A:0.6299 A:0.4807 A:0.4358 A:0.5578 A:0.5185 A:0.5241 A:0.5829
4 30724365 30724365 G/A + rs977931 A:0.1982 A:0.1571 A:0.3204 A:0.1133 A:0.2526 A:0.1723 A:0.1234 A:0.166 A:0.1815 A:0.2013 A:0.1074 A:0.3093 A:0.1365 A:0.126 A:0.1509 A:0.2275
4 40935145 40935145 C/T + rs6846076 T:0.0106 T:0 T:0 T:0 T:0 T:0.009393 T:0.0006285 T:0.001046 T:0.01209 T:0.0005601 T:0 T:0 T:0.001941 T:0.0004066 T:0.0005705 T:5.175e-05
4 70389876 70389876 G/A + rs2988 A:0.7685 A:0.6974 A:0.5942 A:0.6789 A:0.6094 A:0.7562 A:0.6557 A:0.6773 A:0.7667 A:0.7262 A:0.7101 A:0.6253 A:0.736 A:0.6585 A:0.6805 A:0.6232
4 70993963 70993963 G/C + rs9997790 C:0.0015 C:0.0274 C:0 C:0.0249 C:0.0051 C:0.004766 C:0.01802 C:0.01475 C:0.003795 C:0.009619 C:0.02417 C:6.117e-05 C:0.02123 C:0.01997 C:0.01931 C:0.00683
4 78384023 78384023 T/G + rs7677133 G:0.118 G:0.3458 G:0.2173 G:0.2714 G:0.182 G:0.1505 G:0.2741 G:0.2665 G:0.1392 G:0.3581 G:0.2608 G:0.2155 G:0.3114 G:0.2782 G:0.2659 G:0.1848
5 16478091 16478091 G/A + rs162848 A:0.8185 A:0.598 A:0.5526 A:0.7356 A:0.6125 A:0.803 A:0.6984 A:0.653 A:0.8074 A:0.5541 A:0.7153 A:0.5225 A:0.618 A:0.6958 A:0.6766 A:0.6046
5 41004304 41004304 A/T + rs7711722 T:0.5 T:0.4337 T:0.2847 T:0.4185 T:0.2945 T:0.4903 T:0.4367 T:0.41 T:0.5036 T:0.4157 T:0.3878 T:0.2848 T:0.4591 T:0.4376 T:0.4143 T:0.2892
5 73502899 73502899 T/G + rs347233 G:0.0303 G:0.5288 G:0.2619 G:0.3996 G:0.3947 G:0.08401 G:0.4148 G:0.4137 G:0.07169 G:0.6051 G:0.4348 G:0.2669 G:0.3752 G:0.4266 G:0.4177 G:0.4211
5 77349587 77349587 C/T + rs2306342 T:0 T:0.0029 T:0.0298 T:0.001 T:0.0061 T:0.001135 T:0.003605 T:0.004545 T:0.001112 T:0.001937 T:0 T:0.02901 T:0.001121 T:0.003972 T:0.004199 T:0.001463
5 90474182 90474182 G/A + rs2115501 A:0.0825 A:0.0648 A:0.006 A:0.1402 A:0.0481 A:0.08265 A:0.1029 A:0.09021 A:0.07756 A:0.04256 A:0.08536 A:0.002366 A:0.1975 A:0.1109 A:0.08889 A:0.05024
5 96896518 96896518 G/C + rs2548537 C:0.6014 C:0.5793 C:0.5278 C:0.5199 C:0.5276 C:0.5785 C:0.5213 C:0.546 C:0.5735 C:0.6167 C:0.5764 C:0.588 C:0.5245 C:0.5098 C:0.5415 C:0.5718
5 96909814 96909814 A/G + rs2549797 G:0.6021 G:0.5807 G:0.5278 G:0.5199 G:0.5276 G:0.579 G:0.5201 G:0.5469 G:0.5714 G:0.6193 G:0.5996 G:0.5856 G:0.5218 G:0.5077 G:0.5442 G:0.5736
5 136356887 136356887 G/A + rs2546661 A:0.4312 A:0.2075 A:0.3641 A:0.164 A:0.1022 A:0.3681 A:0.1591 A:0.1946 A:0.383 A:0.2397 A:0.2901 A:0.3486 A:0.1789 A:0.1494 A:0.1929 A:0.1106
5 140560648 140560648 G/A + rs250431 A:0.292 A:0.7248 A:0.8254 A:0.7137 A:0.7393 A:0.3706 A:0.699 A:0.6881 A:0.3478 A:0.7184 A:0.653 A:0.8385 A:0.6317 A:0.6936 A:0.6645 A:0.7757
5 141247059 141247059 G/A + rs618096 A:0.7935 A:0.6326 A:0.6597 A:0.5497 A:0.5593 A:0.7544 A:0.5686 A:0.5963 A:0.7591 A:0.6346 A:0.5593 A:0.6198 A:0.5431 A:0.571 A:0.5814 A:0.6055
6 3010225 3010225 C/T + rs2475208 T:0.6256 T:0.4568 T:0.5327 T:0.4871 T:0.5276 T:0.5953 T:0.4773 T:0.4953 T:0.5948 T:0.47 T:0.485 T:0.5082 T:0.5032 T:0.4757 T:0.4933 T:0.5295
6 6151978 6151978 C/G + rs2274394 G:0.2292 G:0.2046 G:0.0923 G:0.2406 G:0.1953 G:0.1993 G:0.232 G:0.2086 G:0.2054 G:0.2033 G:0.2784 G:0.09802 G:0.1817 G:0.2261 G:0.2144 G:0.211
6 26501540 26501540 C/T + rs2295593 T:0.2163 T:0.0605 T:0.122 T:0.1223 T:0.1616 T:0.1877 T:0.1106 T:0.1129 T:0.2104 T:0.06136 T:0.04479 T:0.1098 T:0.08558 T:0.1152 T:0.09758 T:0.16
6 29723936 29723936 A/C + rs2076182 C:0.264 C:0.1643 C:0.2917 C:0.1183 C:0.2883 C:0.214 C:0.1145 C:0.183 C:0.2715 C:0.1993 C:0.1781 C:0.3267 C:0.1114 C:0.1246 C:0.171 C:0.2942
6 29828350 29828350 -/C + rs3215482 C:0.5514 C:0.4726 C:0.6111 C:0.4334 C:0.7526 C:0.5261 C:0.4519 C:0.4969 C:0.5385 C:0.4968 C:0.5678 C:0.6065 C:0.3323 C:0.4454 C:0.4953 C:0.7053
6 30580433 30580433 A/G + rs6902544 G:0.0287 G:0.0086 G:0.003 G:0.0089 G:0.0072 G:0.03517 G:0.009789 G:0.008004 G:0.03114 G:0.01028 G:0.009863 G:0.001806 G:0.0004793 G:0.008213 G:0.007835 G:0.004883
6 31355250 31355250 C/T + rs3819276 T:0.0015 T:0 T:0.0069 T:0.0099 T:0.0286 T:0.004633 T:0.01384 T:0.01254 T:0.004281 T:0.006642 T:0.009556 T:0.006965 T:0.006439 T:0.01422 T:0.0112 T:0.02636
6 31779733 31779733 C/T + rs5030798 T:0.0008 T:0.0058 T:0.0169 T:0.008 T:0.0225 T:0.0009934 T:0.006093 T:0.006892 T:0.0001409 T:0.004213 T:0.03027 T:0.01261 T:0.0005844 T:0.004448 T:0.009249 T:0.01515
6 32155862 32155862 G/A + rs2269425 A:0.0923 A:0.2147 A:0.1706 A:0.1769 A:0.1074 A:0.1052 A:0.1482 A:0.1622 A:0.09716 A:0.3049 A:0.02772 A:0.2048 A:0.2353 A:0.1402 A:0.1414 A:0.08888
6 32660268 32660268 C/G + rs9273472 G:0.649 G:0.8199 G:0.7688 G:0.7266 G:0.8272 G:0.4549 G:0.5251 G:0.5309 G:0.431 G:0.5615 G:0.6198 G:0.6924 G:0.5195 G:0.4449 G:0.5303 G:0.7427
7 12604709 12604709 G/T + rs849778 T:0.7345 T:0.9207 T:0.999 T:0.8996 T:0.9755 T:0.7637 T:0.8906 T:0.9147 T:0.7401 T:0.9315 T:0.8356 T:0.9998 T:0.939 T:0.9004 T:0.8986 T:0.9711
7 20647558 20647558 C/T + rs2893006 T:0.3986 T:0.1671 T:0.0407 T:0.2068 T:0.2117 T:0.3794 T:0.2064 T:0.2064 T:0.3823 T:0.15 T:0.2192 T:0.03408 T:0.2973 T:0.1965 T:0.2041 T:0.245
7 38391829 38391829 A/G + rs6949992 G:0.003 G:0.0663 G:0.1409 G:0.1044 G:0.0542 G:0.01816 G:0.1014 G:0.08498 G:0.01596 G:0.05441 G:0.1588 G:0.1412 G:0.07684 G:0.09783 G:0.09825 G:0.05453
7 45073522 45073522 G/A + rs2289366 A:0.0038 A:0.1138 A:0.005 A:0 A:0.001 A:0.002497 A:0.0005814 A:0.0172 A:0.002228 A:0.1159 A:0 A:0.006672 A:0.000591 A:0.0002159 A:0.01465 A:0.002112
7 73737314 73737314 A/G + rs6460052 G:0.8298 G:0.5389 G:0.5218 G:0.4891 G:0.5511 G:0.768 G:0.4919 G:0.5394 G:0.7752 G:0.5565 G:0.5651 G:0.5265 G:0.4956 G:0.4997 G:0.5454 G:0.5766
7 75418889 75418889 C/T + rs236655 T:0 T:0.0014 T:0.0188 T:0.0159 T:0.0051 T:0.003406 T:0.0156 T:0.01491 T:0.002196 T:0.005263 T:0.008819 T:0.01477 T:0.04009 T:0.0172 T:0.01405 T:0.006886
7 82952543 82952543 C/T + rs976714 T:0.1029 T:0.389 T:0.4405 T:0.3608 T:0.3988 T:0.1379 T:0.3604 T:0.3695 T:0.1299 T:0.4132 T:0.3739 T:0.4938 T:0.391 T:0.3596 T:0.3672 T:0.391
7 87369526 87369526 G/A + rs6966140 A:0.3616 A:0.732 A:0.747 A:0.7763 A:0.7127 A:0.4321 A:0.7726 A:0.7327 A:0.4147 A:0.7113 A:0.7829 A:0.7351 A:0.738 A:0.7753 A:0.7397 A:0.7264
7 87550300 87550300 C/T + rs2032586 T:0.0136 T:0 T:0 T:0 T:0 T:0.008852 T:0 T:0.0008207 T:0.0115 T:0.0005063 T:0 T:0 T:0 T:4.479e-05 T:0.000547 T:3.249e-05
7 100091199 100091199 A/C + rs2272346 C:0.0318 C:0.0418 C:0.0397 C:0.0666 C:0.0552 C:0.03609 C:0.04674 C:0.04188 C:0.03156 C:0.02966 C:0.03089 C:0.04029 C:0.03482 C:0.04655 C:0.04711 C:0.05201
8 3396420 3396420 A/G + rs3802305 G:0.5454 G:0.4957 G:0.251 G:0.5447 G:0.5225 G:0.4692 G:0.4931 G:0.4913 G:0.5208 G:0.4725 G:0.5667 G:0.245 G:0.5347 G:0.5237 G:0.5062 G:0.5104
8 10538544 10538544 A/C + rs6601483 C:0.6399 C:0.6902 C:0.7639 C:0.5964 C:0.6728 C:0.623 C:0.6305 C:0.6486 C:0.6211 C:0.6963 C:0.5199 C:0.761 C:0.6137 C:0.6324 C:0.6231 C:0.6773
8 17114101 17114101 A/C + rs6999470 C:0.0219 C:0.0014 C:0 C:0 C:0 C:0.0134 C:0.0001163 C:0.001271 C:0.01475 C:0.001795 C:0.001941 C:0 C:0 C:2.698e-05 C:0.0007364 C:0
8 22067703 22067703 C/G + rs6997244 G:0.5567 G:0.7824 G:0.9772 G:0.827 G:0.8303 G:0.5912 G:0.8337 G:0.8124 G:0.5856 G:0.7208 G:0.8788 G:0.975 G:0.7873 G:0.8378 G:0.8293 G:0.8349
8 28793781 28793781 A/C + rs2665911 C:0.3011 C:0.5101 C:0.1667 C:0.5537 C:0.3487 C:0.3267 C:0.6299 C:0.4982 C:0.3515 C:0.4721 C:0.5701 C:0.2311 C:0.5566 C:0.5525 C:0.5189 C:0.437
8 65605251 65605251 C/T + rs3765205 T:0.289 T:0.4236 T:0.4355 T:0.3867 T:0.3231 T:0.3034 T:0.3953 T:0.3781 T:0.2974 T:0.4134 T:0.3773 T:0.465 T:0.3187 T:0.3903 T:0.3777 T:0.3312
8 93795937 93795937 A/G + rs3134031 G:0.7625 G:0.67 G:0.7331 G:0.6123 G:0.6564 G:0.7291 G:0.595 G:0.6392 G:0.7337 G:0.74 G:0.6593 G:0.7325 G:0.5737 G:0.5977 G:0.6282 G:0.6235
8 103921800 103921800 G/A + rs6468897 A:0.7511 A:0.4914 A:0.2421 A:0.7097 A:0.5838 A:0.718 A:0.6908 A:0.6337 A:0.7285 A:0.4729 A:0.6578 A:0.236 A:0.7574 A:0.6945 A:0.6567 A:0.6302
8 123534764 123534764 T/C + rs6988591 C:0.2496 C:0.0187 C:0 C:0 C:0 C:0.2233 C:0.0006977 C:0.01577 C:0.22 C:0.01124 C:0 C:0 C:0.0005391 C:0.0005645 C:0.009664 C:0.0004875
8 140514385 140514385 G/A + rs2231522 A:0.0204 A:0 A:0 A:0 A:0 A:0.01634 A:0.0003488 A:0.00113 A:0.01466 A:0.000757 A:0.0004534 A:0 A:0 A:0.0001048 A:0.0006237 A:7.87e-05
9 676954 676954 T/C + rs2296055 C:0.7345 C:0.3689 C:0.5923 C:0.2217 C:0.2802 C:0.6371 C:0.2191 C:0.3085 C:0.6603 C:0.3808 C:0.2578 C:0.5744 C:0.313 C:0.2174 C:0.2951 C:0.2523
9 740901 740901 C/T + rs2296049 T:0.413 T:0.317 T:0.1716 T:0.2117 T:0.184 T:0.3951 T:0.2199 T:0.2427 T:0.409 T:0.3798 T:0.2061 T:0.1938 T:0.1834 T:0.2165 T:0.2479 T:0.1862
9 14759727 14759727 A/G + rs7045154 G:0.9153 G:0.9942 G:1 G:1 G:1 G:0.9294 G:0.9999 G:0.9938 G:0.9188 G:0.9957 G:1 G:1 G:1 G:1 G:0.9966 G:0.9999
9 98209066 98209066 T/C + rs3739666 C:0.7731 C:0.5418 C:0.498 C:0.5706 C:0.5399 C:0.7458 C:0.5605 C:0.5663 C:0.7613 C:0.5006 C:0.5959 C:0.5293 C:0.5752 C:0.5685 C:0.5808 C:0.5353
9 104830901 104830901 A/T + rs4743763 T:0.8427 T:0.4092 T:0.7887 T:0.2614 T:0.2924 T:0.7515 T:0.2797 T:0.3596 T:0.7763 T:0.4057 T:0.3268 T:0.8085 T:0.2058 T:0.2689 T:0.3549 T:0.3024
9 108862979 108862979 G/A + rs3739693 A:0.1793 A:0.0187 A:0.0734 A:0 A:0.0992 A:0.1264 A:0.001279 A:0.02549 A:0.1404 A:0.008841 A:0.01253 A:0.06379 A:0.0002279 A:0.001332 A:0.01181 A:0.07671
9 109037987 109037987 G/A + rs7029890 A:0.211 A:0.0879 A:0 A:0.159 A:0.1196 A:0.1837 A:0.1479 A:0.1323 A:0.1954 A:0.0667 A:0.1506 A:0.001297 A:0.1306 A:0.1554 A:0.1441 A:0.1269
9 110875574 110875574 G/A + rs3739709 A:0.0113 A:0.1643 A:0.1905 A:0.1988 A:0.1524 A:0.04494 A:0.1927 A:0.1774 A:0.03764 A:0.2033 A:0.155 A:0.1825 A:0.1769 A:0.1921 A:0.1754 A:0.1705
9 123373968 123373968 A/G + rs2488599 G:0.7027 G:0.304 G:0.3661 G:0.3221 G:0.4254 G:0.4675 G:0.2469 G:0.3608 G:0.6723 G:0.2885 G:0.4377 G:0.3761 G:0.2747 G:0.3354 G:0.38 G:0.4249
9 128311516 128311516 G/A + rs2231636 A:0.1165 A:0.0072 A:0 A:0 A:0 A:0.06809 A:0.0004651 A:0.005675 A:0.07626 A:0.004348 A:0 A:5.798e-05 A:0 A:0.0004388 A:0.004194 A:0.0003898
X 13662938 13662938 T/C + rs7061978 C:0.1486 C:0.0076 C:0 C:0 C:0 C:0.1356 C:0.0003102 C:0.01129 C:0.1456 C:0.00625 C:0 C:0 C:0 C:0.0003997 C:0.007411 C:0.0004702
X 55488385 55488385 C/T + rs3126255 T:0.0668 T:0.4637 T:0.6217 T:0.6762 T:0.6281 T:0.1304 T:0.6888 T:0.6013 T:0.1103 T:0.4387 T:0.5102 T:0.6349 T:0.7505 T:0.6963 T:0.5777 T:0.6396
X 70352640 70352640 C/T + rs2297871 T:0 T:0 T:0.106 T:0 T:0.0056 T:0 T:0.001189 T:0.008376 T:8.245e-05 T:0.0001921 T:0.004144 T:0.08747 T:0.006 T:0.001435 T:0.004779 T:0.005875
X 91436518 91436518 A/G + rs7050077 G:0.014 G:0.0038 G:0 G:0 G:0 G:0.01278 G:0 G:0.0009344 G:0.01243 G:0.000301 G:0 G:0 G:0 G:7.489e-05 G:0 G:5.224e-05
X 100633567 100633567 G/C + rs2234091 C:0.0279 C:0.0038 C:0 C:0.0013 C:0 C:0.03444 C:0 C:0.002717 C:0.02987 C:0.001986 C:0 C:0 C:0 C:5.859e-05 C:0.001939 C:0
X 124406465 124406465 G/T + rs960869 T:0.3918 T:0.4065 T:0.3573 T:0.4256 T:0.3217 T:0.3974 T:0.4217 T:0.4118 T:0.4038 T:0.4521 T:0.4503 T:0.3425 T:0.441 T:0.4198 T:0.3903 T:0.3362
X 131304109 131304109 G/A + rs7057216 A:0.011 A:0.0019 A:0 A:0 A:0 A:0.01199 A:0 A:0.0007517 A:0.009911 A:0.0003178 A:0 A:0 A:0 A:0 A:0.0005464 A:6.02e-05
X 137030548 137030548 A/G + rs5931046 G:0.2572 G:0.1164 G:0.1348 G:0.171 G:0.3565 G:0.2138 G:0.1476 G:0.168 G:0.2364 G:0.1115 G:0.1941 G:0.1243 G:0.1377 G:0.1522 G:0.1641 G:0.3183
X 152914515 152914515 G/A + rs6526142 A:0.3569 A:0.145 A:0.1296 A:0.3042 A:0.273 A:0.3342 A:0.2929 A:0.2295 A:0.3328 A:0.09453 A:0.1941 A:0.1247 A:0.2558 A:0.2806 A:0.2437 A:0.2396
Y 12720687 12720687 G/T + rs7067496 T:0.931 T:0.1824 T:0.1311 T:0.025 T:0.0692 T:0.69 T:0.05716 T:0.08915 T:0.7451 T:0.09489 T:0.2131 T:0.08762 T:0.002907 T:0.05014 T:0.105 T:0.0332
Y 12739620 12739620 T/C + rs2032599 C:0.0188 C:0.0059 C:0 C:0 C:0 C:0.02452 C:0 C:0.001198 C:0.02299 C:0.0007298 C:0 C:0 C:0 C:6.402e-05 C:0.00142 C:8.686e-05
Y 12776849 12776849 A/C + rs2032600 C:0 C:0.0412 C:0 C:0.0042 C:0 C:0.001751 C:0.004274 C:0.004949 C:0.00131 C:0.02136 C:0.004317 C:0 C:0.0001729 C:0.003886 C:0.03123 C:8.61e-05
Y 12786229 12786229 G/A + rs20320 A:0 A:0.0235 A:0 A:0.0167 A:0 A:0.001751 A:0.006944 A:0.008052 A:0.0003285 A:0.01438 A:0.03539 A:0 A:0.0001714 A:0.007779 A:0.01915 A:0.006652
Y 12790481 12790481 G/A + rs20321 A:0.0031 A:0 A:0 A:0 A:0 A:0.003503 A:0.03579 A:0.009614 A:0.006544 A:0.0008574 A:0.000392 A:0 A:0.0001708 A:0.01933 A:0.007725 A:8.598e-05
Y 13357844 13357844 G/A + rs9341278 A:0 A:0 A:0.0328 A:0.0958 A:0.0038 A:0 A:0.005876 A:0.07612 A:0.0003937 A:0.001661 A:0.0004146 A:0.0409 A:0.6944 A:0.01712 A:0.07514 A:0.00051
Y 13479657 13479657 G/C + rs2032653 C:0.8966 C:0.1765 C:0.082 C:0.025 C:0 C:0.6567 C:0.05662 C:0.07925 C:0.7276 C:0.09543 C:0.2131 C:0.01325 C:0.002578 C:0.04961 C:0.1046 C:0.01038
Y 19707378 19707378 G/A + rs2032620 A:0.0094 A:0 A:0 A:0 A:0 A:0.02977 A:0 A:0.001341 A:0.02715 A:0.0007612 A:0 A:0.0004965 A:0 A:3.679e-05 A:0.001555 A:0
Y 19731995 19731995 A/C + rs2032672 C:0 C:0.0235 C:0 C:0.0167 C:0 C:0.001751 C:0.006944 C:0.008169 C:0.0003272 C:0.01461 C:0.03532 C:0 C:0.0001713 C:0.008068 C:0.01758 C:0.006622
};
  }
  elsif($assembly eq 'GRCh37') {
    print OUT qq{1 1182031 1182031 C/T + rs7539412 T:0.0439 T:0.0101 T:0 T:0 T:0 T:0.0364 T:0.0007292 T:0.001832 T:0.1 T:0 T:0 T:0 T:0 T:0 T:0 T:0
1 1686962 1686962 C/T + rs2076327 T:0.1846 T:0.4078 T:0.3145 T:0.493 T:0.3804 T:0.2238 T:0.4971 T:0.4491 T:0.2125 T:0.4256 T:0.5762 T:0.2899 T:0.492 T:0.5046 T:0.4851 T:0.4091
1 8390495 8390495 A/G + rs1058619 G:0.7368 G:0.7824 G:0.8482 G:0.7366 G:0.7924 G:0.7468 G:0.742 G:0.7694 G:0.7563 G:0.7774 G:0.6359 G:0.8537 G:0.8548 G:0.7484 G:0.7389 G:0.7807
1 11105545 11105545 T/C + rs2273343 C:0 C:0 C:0.0268 C:0 C:0.001 C:0 C:0.0002326 C:0.002332 C:0.0001999 C:2.984e-05 C:0 C:0.02913 C:0 C:0.0001087 C:0.002025 C:0.0013
1 11114940 11114940 G/T + rs13616 T:0.2481 T:0.1167 T:0.1806 T:0.0249 T:0.229 T:0.2072 T:0.02267 T:0.08616 T:0.2176 T:0.1433 T:0.05073 T:0.1767 T:0.05429 T:0.02378 T:0.07749 T:0.1668
1 11796240 11796240 C/T + rs6668392 T:0 T:0.0216 T:0 T:0.0338 T:0.0092 T:0.003876 T:0.03116 T:0.02254 T:0.004753 T:0.007478 T:0.0201 T:0 T:0.04682 T:0.0325 T:0.01842 T:0.009781
1 12918819 12918819 C/T + rs9659511 T:0.1346 T:0.1124 T:0.3214 T:0.1421 T:0.1442 T:0.142 T:0.1364 T:0.1533 T:0.1474 T:0.1043 T:0.1157 T:0.3109 T:0.1857 T:0.146 T:0.1434 T:0.1286
1 18808292 18808292 C/A + rs2992753 A:0.6604 A:0.7277 A:0.8958 A:0.6113 A:0.4847 A:0.665 A:0.6287 A:0.6638 A:0.67 A:0.8076 A:0.6959 A:0.8954 A:0.6439 A:0.6284 A:0.6655 A:0.5068
1 19487656 19487656 T/C + rs3748759 C:0.5197 C:0.5548 C:0.7986 C:0.4294 C:0.5184 C:0.5291 C:0.43 C:0.4896 C:0.5252 C:0.5752 C:0.3456 C:0.8051 C:0.3984 C:0.4438 C:0.4506 C:0.4848
1 19582563 19582563 G/A + rs1557133 A:0.9206 A:0.8833 A:0.9266 A:0.8549 A:0.8978 A:0.9199 A:0.8373 A:0.8635 A:0.9143 A:0.8835 A:0.7688 A:0.9166 A:0.8747 A:0.8479 A:0.8523 A:0.8674
10 3186458 3186458 C/G + rs3752803 G:0.6218 G:0.549 G:0.6696 G:0.508 G:0.5368 G:0.6291 G:0.5068 G:0.5328 G:0.636 G:0.6149 G:0.4908 G:0.6717 G:0.4657 G:0.49 G:0.5214 G:0.53
10 7785087 7785087 G/A + rs3736966 A:0.0212 A:0.0375 A:0.2411 A:0.0626 A:0.1084 A:0.03065 A:0.06446 A:0.08329 A:0.02488 A:0.03598 A:0.09121 A:0.2536 A:0.09442 A:0.06703 A:0.08262 A:0.1145
10 13340293 13340293 A/G + rs7916926 G:0.0499 G:0.5648 G:0.6786 G:0.4215 G:0.4836 G:0.1243 G:0.4337 G:0.4639 G:0.1035 G:0.5885 G:0.3889 G:0.6955 G:0.4541 G:0.4454 G:0.4627 G:0.4695
10 15573054 15573054 T/C + rs9333241 C:0.0242 C:0.0029 C:0 C:0 C:0 C:0.02111 C:0 C:0.001647 C:0.02284 C:0.001238 C:0 C:0 C:0 C:2.913e-05 C:0.0008569 C:0
10 16797033 16797033 A/G + rs7910261 G:0.8911 G:0.3127 G:0.4048 G:0.2654 G:0.2894 G:0.7826 G:0.2526 G:0.295 G:0.81 G:0.2237 G:0.2486 G:0.415 G:0.2222 G:0.255 G:0.289 G:0.2539
10 25140411 25140411 G/C + rs1326192 C:0.3374 C:0.353 C:0.506 C:0.4722 C:0.545 C:0.3593 C:0.4688 C:0.4552 C:0.3637 C:0.3617 C:0.5291 C:0.5229 C:0.3622 C:0.4686 C:0.4711 C:0.5452
10 28101455 28101455 C/A + rs3737184 A:0.0015 A:0.0908 A:0.1776 A:0.0089 A:0.1329 A:0.003177 A:0.00186 A:0.05041 A:0.002358 A:0.1118 A:0.002339 A:0.2236 A:0.03742 A:0.002651 A:0.03611 A:0.1126
10 43316170 43316170 A/G + rs3982211 G:0.6271 G:0.4885 G:0.626 G:0.4185 G:0.589 G:0.5618 G:0.3788 G:0.4599 G:0.5939 G:0.5122 G:0.4353 G:0.6343 G:0.3903 G:0.3888 G:0.4394 G:0.5597
10 52576068 52576068 G/A + rs4245008 A:0.5946 A:0.889 A:0.9663 A:0.838 A:0.7791 A:0.6312 A:0.8508 A:0.8336 A:0.6184 A:0.9128 A:0.7954 A:0.9808 A:0.7715 A:0.8524 A:0.8303 A:0.7623
10 75095253 75095253 G/A + rs7910541 A:0.177 A:0.0994 A:0.3105 A:0.0547 A:0.226 A:0.1496 A:0.06663 A:0.1096 A:0.1586 A:0.08698 A:0.1167 A:0.3029 A:0.04417 A:0.06152 A:0.08687 A:0.2263
11 2339195 2339195 G/A + rs2234316 A:0.0333 A:0 A:0 A:0 A:0 A:0.03247 A:0.0001163 A:0.00204 A:0.02874 A:0.0007769 A:0 A:6.259e-05 A:0 A:0.0001309 A:0.00141 A:0
11 4615234 4615234 A/G + rs4910511 G:0.0265 G:0.245 G:0.4405 G:0.34 G:0.3487 G:0.08178 G:0.3285 G:0.3242 G:0.06484 G:0.3198 G:0.2778 G:0.4075 G:0.3542 G:0.3365 G:0.3097 G:0.3626
11 4615987 4615987 A/C + rs1876025 C:0.0537 C:0.0043 C:0 C:0.001 C:0 C:0.04634 C:0.0009307 C:0.003749 C:0.0524 C:0.001846 C:0.0003046 C:0 C:0 C:0.0004031 C:0.001458 C:9.746e-05
11 5068177 5068177 A/T + rs2500018 T:0.5915 T:0.4971 T:0.2619 T:0.508 T:0.5521 T:0.5552 T:0.5007 T:0.4925 T:0.5654 T:0.4971 T:0.5515 T:0.2763 T:0.4809 T:0.5001 T:0.5062 T:0.5317
11 7063755 7063755 C/T + rs7123944 T:0.8533 T:0.9885 T:1 T:0.998 T:1 T:0.8719 T:0.9977 T:0.9885 T:0.8613 T:0.9896 T:0.985 T:1 T:1 T:0.9987 T:0.9896 T:0.9998
11 7716855 7716855 T/A + rs4528317 A:0.8517 A:0.9697 A:1 A:0.998 A:1 A:0.8723 A:0.9993 A:0.9893 A:0.8609 A:0.991 A:0.9962 A:1 A:1 A:0.999 A:0.9901 A:0.9996
11 7817538 7817538 T/C + rs7949771 C:0.3472 C:0.2911 C:0.128 C:0.2763 C:0.1196 C:0.3231 C:0.2819 C:0.2557 C:0.3451 C:0.2861 C:0.2303 C:0.1323 C:0.2788 C:0.2801 C:0.2504 C:0.151
11 8118200 8118200 C/A + rs7931842 A:0.9856 A:0.9986 A:1 A:1 A:1 A:0.9939 A:1 A:0.9994 A:0.992 A:0.9996 A:1 A:1 A:1 A:1 A:0.9993 A:1
11 8284860 8284860 G/T + rs204930 T:0.1339 T:0.255 T:0.2351 T:0.2376 T:0.1769 T:0.1324 T:0.2239 T:0.2225 T:0.1472 T:0.2294 T:0.263 T:0.227 T:0.2655 T:0.2309 T:0.2376 T:0.1775
11 9096086 9096086 G/A + rs2003885 A:0.1051 A:0.0072 A:0 A:0 A:0 A:0.092 A:0.0002328 A:0.006094 A:0.08758 A:0.003722 A:0 A:0 A:0 A:0.0001523 A:0.002371 A:0.0001624
12 6030301 6030301 G/A + rs3741901 A:0.0045 A:0.0389 A:0.0139 A:0.1064 A:0.1166 A:0.02044 A:0.12 A:0.0917 A:0.01707 A:0.05393 A:0.1156 A:0.02442 A:0.07702 A:0.1136 A:0.09967 A:0.1299
12 6882121 6882121 G/- + rs3214312 -:0.0061 -:0.0461 -:0.1359 -:0.1044 -:0.0215 -:0.02205 -:0.09886 -:0.08039 -:0.01637 -:0.05027 -:0.05388 -:0.1252 -:0.1257 -:0.1003 -:0.08824 -:0.02309
12 12483275 12483275 G/A + rs1861676 A:0.5454 A:0.6628 A:0.38 A:0.7266 A:0.7117 A:0.581 A:0.7542 A:0.6869 A:0.5701 A:0.6349 A:0.6675 A:0.3754 A:0.7124 A:0.7475 A:0.7013 A:0.7418
12 21331987 21331987 C/T + rs2291076 T:0.0477 T:0.4035 T:0.245 T:0.4642 T:0.4223 T:0.1089 T:0.4601 T:0.4144 T:0.0869 T:0.4531 T:0.4287 T:0.2541 T:0.4096 T:0.4555 T:0.4094 T:0.4765
12 22814036 22814036 A/G + rs4963793 G:0.0794 G:0.196 G:0.6925 G:0.1292 G:0.2178 G:0.0892 G:0.1159 G:0.1969 G:0.08337 G:0.2878 G:0.1307 G:0.6998 G:0.2134 G:0.1118 G:0.1838 G:0.1949
12 25147300 25147300 G/C + rs699030 C:0.8593 C:0.9914 C:1 C:1 C:1 C:0.8847 C:0.9999 C:0.9911 C:0.8743 C:0.9942 C:1 C:1 C:1 C:0.9997 C:0.9946 C:0.9998
12 31240960 31240960 G/A + rs7309189 A:0.7088 A:0.6066 A:0.8095 A:0.4314 A:0.5746 A:0.685 A:0.451 A:0.52 A:0.6832 A:0.6578 A:0.394 A:0.8083 A:0.4148 A:0.4338 A:0.4887 A:0.5643
12 32490764 32490764 G/A + rs2270786 A:0.0242 A:0.1427 A:0.3036 A:0.2276 A:0.0941 A:0.04975 A:0.207 A:0.1841 A:0.04361 A:0.1257 A:0.161 A:0.2886 A:0.2881 A:0.2141 A:0.1902 A:0.08208
12 39233771 39233771 A/G + rs3803020 G:0.0711 G:0.5 G:0.1855 G:0.495 G:0.4049 G:0.1281 G:0.4987 G:0.4306 G:0.1185 G:0.4506 G:0.4507 G:0.1491 G:0.4581 G:0.5046 G:0.4501 G:0.4224
12 40499347 40499347 G/A + rs7305377 A:0.0212 A:0.0908 A:0.0308 A:0.1581 A:0.0511 A:0.03866 A:0.1321 A:0.1281 A:0.0366 A:0.1136 A:0.08116 A:0.02681 A:0.2627 A:0.1481 A:0.1442 A:0.08264
13 21215379 21215379 A/G + rs7332496 G:0.0098 G:0 G:0 G:0 G:0 G:0.007431 G:0 G:0.0007988 G:0.01006 G:0.0007105 G:0 G:0 G:0 G:2.953e-05 G:0.0008251 G:0
13 25451375 25451375 A/G + rs3803215 G:0.0469 G:0.0274 G:0.0764 G:0.002 G:0.0031 G:0.04607 G:0.001047 G:0.01437 G:0.04088 G:0.03098 G:0.001279 G:0.08681 G:0.004865 G:0.001327 G:0.00918 G:0.003106
13 27829310 27829310 G/A + rs3094293 A:0.944 A:0.9914 A:1 A:1 A:1 A:0.9561 A:0.9995 A:0.996 A:0.9443 A:0.9977 A:0.9978 A:1 A:1 A:0.9998 A:0.9974 A:0.9999
13 36744910 36744910 C/T + rs2296968 T:0.1324 T:0.3559 T:0.4911 T:0.1769 T:0.227 T:0.1505 T:0.1969 T:0.2425 T:0.1426 T:0.3799 T:0.134 T:0.4954 T:0.1642 T:0.2013 T:0.2176 T:0.2452
13 45914957 45914957 T/C + rs2234216 C:0.9085 C:0.5749 C:0.3353 C:0.6223 C:0.4632 C:0.8484 C:0.6346 C:0.5878 C:0.8622 C:0.4781 C:0.6003 C:0.3306 C:0.6025 C:0.632 C:0.5657 C:0.487
13 76055397 76055397 C/T + rs9573565 T:0.1165 T:0.2594 T:0.2103 T:0.2922 T:0.1237 T:0.1406 T:0.2672 T:0.2266 T:0.1507 T:0.1778 T:0.3385 T:0.1852 T:0.1881 T:0.2763 T:0.2635 T:0.1505
13 111164328 111164328 G/A + rs7320105 A:0.0053 A:0.0014 A:0.001 A:0.002 A:0 A:0.004771 A:0.002958 A:0.001904 A:0.004912 A:0.0005064 A:0.004679 A:0 A:0.001211 A:0.002611 A:0.002012 A:6.498e-05
13 113460624 113460624 C/T + rs2275282 T:0 T:0 T:0.0268 T:0 T:0.002 T:0.0004539 T:0 T:0.001359 T:0.0002627 T:9.35e-05 T:0 T:0.01559 T:0 T:2.705e-05 T:0.002039 T:0.001421
13 113819009 113819009 T/A + rs494860 A:0.0605 A:0.2507 A:0.3313 A:0.2256 A:0.5143 A:0.08012 A:0.2283 A:0.2646 A:0.07693 A:0.2773 A:0.2043 A:0.3588 A:0.2455 A:0.2287 A:0.2451 A:0.4521
14 21945415 21945415 A/G + rs2282046 G:0.4531 G:0.4236 G:0.2649 G:0.4364 G:0.4223 G:0.424 G:0.4372 G:0.4203 G:0.4305 G:0.4758 G:0.6192 G:0.2538 G:0.3566 G:0.4256 G:0.4485 G:0.4061
14 21968876 21968876 G/A + rs3752411 A:0.146 A:0.2493 A:0.1042 A:0.1471 A:0.1186 A:0.1311 A:0.1388 A:0.15 A:0.1362 A:0.3206 A:0.09773 A:0.1032 A:0.1033 A:0.134 A:0.1477 A:0.1045
14 22409571 22409571 A/G + rs2178779 G:0.6649 G:0.5159 G:0.3383 G:0.4851 G:0.4499 G:0.6533 G:0.5078 G:0.5001 G:0.659 G:0.5292 G:0.4913 G:0.3353 G:0.5291 G:0.508 G:0.5051 G:0.4494
14 23549785 23549785 T/C + rs3811182 C:0.6914 C:0.3689 C:0.4554 C:0.3777 C:0.5102 C:0.6512 C:0.412 C:0.4323 C:0.6601 C:0.3541 C:0.4261 C:0.4271 C:0.3503 C:0.4166 C:0.4112 C:0.529
14 23829253 23829253 G/A + rs2231803 A:0.0348 A:0.0029 A:0 A:0 A:0 A:0.02542 A:0.0003488 A:0.001804 A:0.02488 A:0.001104 A:0 A:5.8e-05 A:0 A:0.0001442 A:0.0007326 A:0.0001304
14 24906515 24906515 C/T + rs3742521 T:0.3933 T:0.1138 T:0.125 T:0.0527 T:0.0665 T:0.3311 T:0.0407 T:0.07816 T:0.3589 T:0.07701 T:0.04504 T:0.1299 T:0.03799 T:0.04396 T:0.06504 T:0.07048
14 52509483 52509483 T/C + rs4898729 C:0.9039 C:0.7867 C:0.8542 C:0.7873 C:0.8998 C:0.8724 C:0.7669 C:0.805 C:0.8784 C:0.8035 C:0.7561 C:0.841 C:0.8329 C:0.7654 C:0.8091 C:0.8877
14 55864130 55864130 A/G + rs8003279 G:0.1142 G:0.1902 G:0.1617 G:0.329 G:0.274 G:0.1492 G:0.3212 G:0.2766 G:0.1452 G:0.1601 G:0.1904 G:0.1822 G:0.4481 G:0.3153 G:0.276 G:0.2864
14 64421392 64421392 G/A + rs1255874 A:0.6157 A:0.7291 A:0.5982 A:0.6899 A:0.7699 A:0.6174 A:0.7016 A:0.7048 A:0.618 A:0.7607 A:0.7054 A:0.6059 A:0.7123 A:0.6933 A:0.7061 A:0.7759
14 92482254 92482254 T/C + rs7158303 C:0.2133 C:0.0159 C:0.0169 C:0.001 C:0.0491 C:0.144 C:0.0009352 C:0.02158 C:0.1693 C:0.009539 C:0.0003463 C:0.0258 C:0 C:0.001228 C:0.01824 C:0.04351
15 34103272 34103272 C/G + rs4780181 G:0.8222 G:0.8501 G:0.9484 G:0.8101 G:0.8282 G:0.8405 G:0.7864 G:0.8316 G:0.8346 G:0.8681 G:0.7293 G:0.9561 G:0.883 G:0.7952 G:0.8354 G:0.8405
15 39880941 39880941 G/T + rs2292304 T:0 T:0 T:0.1796 T:0.001 T:0.0215 T:0.0002273 T:0.0005818 T:0.02003 T:0.0007825 T:0.0003708 T:0.0007758 T:0.1982 T:0.01888 T:0.0006439 T:0.01234 T:0.01578
15 41021051 41021051 T/A + rs7161941 A:0.1082 A:0.0058 A:0 A:0 A:0 A:0.08557 A:0 A:0.006252 A:0.09243 A:0.002591 A:0 A:0 A:0 A:0.0001791 A:0.002735 A:9.747e-05
15 41146581 41146581 T/C + rs3736287 C:0.6573 C:0.8184 C:0.7143 C:0.7863 C:0.771 C:0.6882 C:0.7781 C:0.7752 C:0.6716 C:0.8452 C:0.7385 C:0.6921 C:0.8313 C:0.7805 C:0.7685 C:0.7693
15 51204395 51204395 -/T + rs3214932 T:0.6384 T:0.611 T:0.375 T:0.5338 T:0.6544 T:0.6231 T:0.5419 T:0.5721 T:0.6244 T:0.6729 T:0.6512 T:0.4287 T:0.5703 T:0.5398 T:0.5734 T:0.609
15 51634255 51634255 T/G + rs2446421 G:0.9115 G:0.6254 G:0.5476 G:0.4811 G:0.4949 G:0.8292 G:0.4549 G:0.5016 G:0.8506 G:0.6137 G:0.4217 G:0.5831 G:0.3243 G:0.4413 G:0.4787 G:0.4913
15 54841874 54841874 A/G + rs9302181 G:0.1346 G:0.2853 G:0.5258 G:0.3887 G:0.5368 G:0.1775 G:0.3804 G:0.3724 G:0.1646 G:0.2679 G:0.3257 G:0.5288 G:0.4178 G:0.3851 G:0.358 G:0.4436
15 55669169 55669169 A/G + rs2289325 G:0 G:0 G:0.0099 G:0.0089 G:0.002 G:0 G:0.001104 G:0.006159 G:0 G:0 G:0.0008128 G:0.01322 G:0.04857 G:0.001604 G:0.004941 G:0.0002599
15 59963561 59963561 C/T + rs6151546 T:0.0197 T:0.0029 T:0 T:0 T:0 T:0.01837 T:0 T:0.001462 T:0.02265 T:0.0006911 T:0 T:0 T:0 T:0 T:0.0007841 T:0
15 62244026 62244026 T/C + rs8026956 C:0.034 C:0.0014 C:0 C:0 C:0 C:0.02838 C:0 C:0.002103 C:0.0317 C:0.0003688 C:0 C:0 C:0 C:5.463e-05 C:0.001117 C:0
16 1573810 1573810 A/G + rs2745176 G:0.9985 G:0.9827 G:1 G:0.9503 G:0.9857 G:0.9902 G:0.9505 G:0.9644 G:0.9924 G:0.9836 G:0.9497 G:0.9999 G:0.9382 G:0.9525 G:0.9588 G:0.9755
16 1962046 1962046 C/T + rs1742399 T:0.2345 T:0.3386 T:0.4931 T:0.3837 T:0.589 T:0.22 T:0.3418 T:0.3903 T:0.2372 T:0.3021 T:0.328 T:0.4985 T:0.4916 T:0.3687 T:0.3936 T:0.5626
16 2549969 2549969 C/T + rs9938971 T:0.0756 T:0.0043 T:0 T:0 T:0 T:0.05937 T:0.0001192 T:0.004821 T:0.06657 T:0.003634 T:0.0001016 T:0 T:0 T:0.000189 T:0.002924 T:0.0002274
16 3065924 3065924 G/A + rs2269911 A:0.1399 A:0.1585 A:0.4028 A:0.1978 A:0.2597 A:0.1385 A:0.1903 A:0.2188 A:0.1468 A:0.1865 A:0.1564 A:0.3826 A:0.2588 A:0.1988 A:0.2068 A:0.2632
16 4833615 4833615 C/T + rs6500634 T:0.0885 T:0.1657 T:0.0694 T:0.2197 T:0.2924 T:0.1022 T:0.1671 T:0.1645 T:0.09582 T:0.09632 T:0.2729 T:0.07353 T:0.1243 T:0.175 T:0.1731 T:0.2787
16 11773743 11773743 G/T + rs7188511 T:0.4342 T:0.3285 T:0.4286 T:0.5547 T:0.4724 T:0.44 T:0.5441 T:0.4954 T:0.4404 T:0.2977 T:0.704 T:0.4137 T:0.5338 T:0.5564 T:0.5141 T:0.4987
16 19644524 19644524 T/C + rs8056927 C:0.2352 C:0.0187 C:0 C:0.002 C:0 C:0.1607 C:0.003023 C:0.01487 C:0.1774 C:0.01182 C:0.006095 C:0 C:0.003012 C:0.00293 C:0.01496 C:0.0004225
16 20043330 20043330 T/C + rs2147865 C:0.6029 C:0.4092 C:0.7827 C:0.328 C:0.4029 C:0.5447 C:0.3034 C:0.3913 C:0.557 C:0.3694 C:0.3057 C:0.7882 C:0.4211 C:0.3126 C:0.3693 C:0.4051
16 23999783 23999783 C/T + rs3729887 T:0.0144 T:0.0187 T:0 T:0.0209 T:0.0041 T:0.01775 T:0.01919 T:0.01866 T:0.0143 T:0.009432 T:0.01041 T:5.849e-05 T:0.05669 T:0.02041 T:0.01887 T:0.01015
16 28618446 28618446 T/C + rs2925623 C:0.23 C:0.4784 C:0.2639 C:0.335 C:0.2239 C:0.2003 C:0.3258 C:0.3572 C:0.232 C:0.4758 C:0.3006 C:0.2024 C:0.4608 C:0.3787 C:0.3641 C:0.2346
17 1371088 1371088 T/C + rs8069059 C:0.9244 C:0.7075 C:0.5526 C:0.7018 C:0.729 C:0.8871 C:0.7336 C:0.7464 C:0.901 C:0.7553 C:0.7582 C:0.6524 C:0.6847 C:0.7529 C:0.7558 C:0.7452
17 2186100 2186100 C/T + rs749240 T:0.6914 T:0.3573 T:0.2371 T:0.3618 T:0.3783 T:0.6096 T:0.3493 T:0.3635 T:0.6356 T:0.335 T:0.3082 T:0.2644 T:0.3361 T:0.3576 T:0.3569 T:0.3744
17 3962450 3962450 A/G + rs9900661 G:0.292 G:0.1657 G:0.1835 G:0.1849 G:0.1779 G:0.2749 G:0.1833 G:0.1856 G:0.2818 G:0.173 G:0.1435 G:0.1707 G:0.1971 G:0.1792 G:0.1946 G:0.1745
17 6483046 6483046 T/C + rs3744720 C:0.8253 C:0.6801 C:0.7093 C:0.5467 C:0.637 C:0.7821 C:0.5747 C:0.6218 C:0.7833 C:0.7462 C:0.5103 C:0.7003 C:0.5492 C:0.5698 C:0.6166 C:0.6376
17 6913718 6913718 C/T + rs2229172 T:0.0045 T:0 T:0 T:0 T:0 T:0.0009079 T:0 T:0.000199 T:0.00281 T:8.935e-05 T:0 T:0 T:0 T:1.791e-05 T:0.0001823 T:0
17 11584124 11584124 A/G + rs9916482 G:0.053 G:0.0029 G:0 G:0 G:0 G:0.05175 G:0.0004651 G:0.004383 G:0.05667 G:0.003403 G:0.006101 G:0 G:0 G:0.0001977 G:0.002195 G:0.0001303
17 11648444 11648444 T/C + rs3744580 C:0.4228 C:0.7565 C:0.6706 C:0.6918 C:0.6902 C:0.4505 C:0.6931 C:0.6879 C:0.4436 C:0.8174 C:0.685 C:0.6771 C:0.6664 C:0.6932 C:0.7057 C:0.6728
17 15848637 15848637 C/T + rs2228101 T:0.0983 T:0.0058 T:0 T:0.002 T:0 T:0.08023 T:0.0003491 T:0.005913 T:0.08377 T:0.004136 T:0 T:6.399e-05 T:0 T:0.0004702 T:0.004254 T:7.021e-05
17 26962640 26962640 G/A + rs665835 A:0.0008 A:0.0735 A:0.0327 A:0 A:0.002 A:0.001365 A:0.0008141 A:0.01778 A:0.001934 A:0.1068 A:0 A:0.03732 A:0 A:0.0001859 A:0.01348 A:0.0007638
17 29486152 29486152 G/A + rs2952976 A:0.23 A:0.5288 A:0.5486 A:0.7157 A:0.6227 A:0.31 A:0.7095 A:0.6268 A:0.2908 A:0.5061 A:0.7772 A:0.5606 A:0.6599 A:0.6971 A:0.6677 A:0.6303
18 9775372 9775372 A/G + rs6506697 G:0.8487 G:0.7118 G:0.6647 G:0.6173 G:0.5961 G:0.8125 G:0.6216 G:0.6642 G:0.8209 G:0.7738 G:0.605 G:0.6704 G:0.6638 G:0.624 G:0.6758 G:0.6256
18 13072979 13072979 T/C + rs1787013 C:0.5 C:0.4553 C:0.3571 C:0.4583 C:0.4162 C:0.5014 C:0.4494 C:0.4508 C:0.4989 C:0.4713 C:0.4938 C:0.3501 C:0.4522 C:0.4609 C:0.4537 C:0.4084
18 30804756 30804756 C/T + rs466113 T:0.9826 T:1 T:1 T:1 T:1 T:0.9859 T:1 T:0.9991 T:0.9866 T:0.9996 T:1 T:1 T:1 T:1 T:1 T:1
18 33557466 33557466 A/G + rs2276314 G:0.323 G:0.1859 G:0.248 G:0.2137 G:0.2352 G:0.3155 G:0.2229 G:0.2286 G:0.3197 G:0.1736 G:0.2857 G:0.253 G:0.2216 G:0.2231 G:0.2159 G:0.2382
18 34340651 34340651 C/T + rs9946890 T:0.1808 T:0.0144 T:0.001 T:0.003 T:0 T:0.15 T:0.0005815 T:0.0114 T:0.1614 T:0.006829 T:0.001427 T:0.0001747 T:0 T:0.0005439 T:0.007896 T:0.000195
18 43206985 43206985 A/G + rs1484873 G:0.8714 G:0.8818 G:0.8532 G:0.9314 G:0.9315 G:0.8833 G:0.9484 G:0.91 G:0.8765 G:0.832 G:0.973 G:0.8625 G:0.8426 G:0.9453 G:0.9203 G:0.9376
18 47787402 47787402 G/A + rs1899671 A:0.4788 A:0.7421 A:0.7917 A:0.7127 A:0.7076 A:0.5122 A:0.692 A:0.699 A:0.5075 A:0.7604 A:0.7143 A:0.797 A:0.7096 A:0.6883 A:0.7091 A:0.7195
18 61264298 61264298 G/A + rs1020694 A:0.9554 A:0.7925 A:0.875 A:0.8072 A:0.819 A:0.9301 A:0.7917 A:0.8031 A:0.9373 A:0.7671 A:0.6934 A:0.8784 A:0.8236 A:0.796 A:0.7767 A:0.7844
18 77246978 77246978 G/A + rs754096 A:0.1891 A:0.4438 A:0.4653 A:0.507 A:0.5072 A:0.2178 A:0.4773 A:0.4779 A:0.2218 A:0.4279 A:0.602 A:0.4674 A:0.5152 A:0.5087 A:0.5207 A:0.5184
19 1062164 1062164 T/C + rs4147920 C:0.1747 C:0.0202 C:0.1736 C:0.0408 C:0.0583 C:0.1486 C:0.03938 C:0.05761 C:0.1532 C:0.02024 C:0.02796 C:0.1638 C:0.056 C:0.0412 C:0.04161 C:0.06411
19 1123652 1123652 G/A + rs4807570 A:0.289 A:0.3761 A:0.1141 A:0.1998 A:0.2362 A:0.2794 A:0.1992 A:0.2464 A:0.3107 A:0.377 A:0.2039 A:0.1345 A:0.2528 A:0.2228 A:0.2287 A:0.2296
19 4816241 4816241 G/A + rs1046673 A:0.2935 A:0.1311 A:0.0675 A:0.1501 A:0.1258 A:0.2851 A:0.1555 A:0.1432 A:0.2906 A:0.09902 A:0.1821 A:0.08516 A:0.07281 A:0.1506 A:0.1502 A:0.1281
19 7687446 7687446 G/A + rs4134852 A:0.0424 A:0.0029 A:0.001 A:0 A:0 A:0.04607 A:0 A:0.003292 A:0.04834 A:0.001727 A:0 A:5.802e-05 A:0 A:8.97e-06 A:0.001094 A:0.0001624
19 8654916 8654916 G/A + rs3923268 A:0.0045 A:0.2277 A:0.4802 A:0.0964 A:0.2403 A:0.01863 A:0.07419 A:0.1525 A:0.01541 A:0.2389 A:0.1627 A:0.5215 A:0.08333 A:0.07582 A:0.1322 A:0.2283
19 9076929 9076929 T/G + rs2547075 G:0.2005 G:0.2118 G:0.249 G:0.3121 G:0.272 G:0.2119 G:0.2883 G:0.2628 G:0.2074 G:0.2172 G:0.1747 G:0.2695 G:0.24 G:0.2905 G:0.2529 G:0.2824
19 10668673 10668673 A/G + rs1982074 G:0.0378 G:0.33 G:0.2192 G:0.163 G:0.227 G:0.06083 G:0.194 G:0.2059 G:0.05321 G:0.3352 G:0.2982 G:0.211 G:0.1475 G:0.1795 G:0.213 G:0.2451
19 11350433 11350433 T/C + rs1541922 C:0.1112 C:0.0072 C:0.0079 C:0.001 C:0.001 C:0.08588 C:0.0003623 C:0.007106 C:0.09232 C:0.006168 C:0.0001028 C:0.00382 C:0 C:0.0003659 C:0.004264 C:0.0006919
19 13029188 13029188 G/A + rs2965214 A:0.5825 A:0.5504 A:0.8482 A:0.6869 A:0.7301 A:0.5932 A:0.654 A:0.6585 A:0.587 A:0.5502 A:0.724 A:0.8499 A:0.7031 A:0.6383 A:0.6611 A:0.7233
19 13264398 13264398 C/T + rs1042164 T:0.0159 T:0.1124 T:0.0129 T:0.175 T:0.0654 T:0.04172 T:0.1673 T:0.1346 T:0.0354 T:0.08011 T:0.06972 T:0.01963 T:0.301 T:0.1744 T:0.1269 T:0.07583
2 27424603 27424603 C/T + rs1064845 T:0.0242 T:0.0014 T:0 T:0 T:0 T:0.02088 T:0 T:0.001539 T:0.02326 T:0.000536 T:0 T:0 T:0 T:8.953e-06 T:0.0005468 T:3.249e-05
2 27481652 27481652 A/C + rs1992290 C:0.4183 C:0.0519 C:0 C:0.0129 C:0.002 C:0.3706 C:0.01628 C:0.03649 C:0.4003 C:0.02201 C:0.03902 C:0 C:0.008388 C:0.01201 C:0.02445 C:0.003249
2 37456032 37456032 C/T + rs2098386 T:0.9902 T:0.9496 T:1 T:0.8549 T:0.8906 T:0.9648 T:0.8358 T:0.8805 T:0.9727 T:0.9504 T:0.9244 T:0.9995 T:0.8408 T:0.8382 T:0.8863 T:0.8587
2 43903483 43903483 T/C + rs4952995 C:0.6142 C:0.5634 C:0.2619 C:0.7187 C:0.6309 C:0.6489 C:0.7245 C:0.6355 C:0.6219 C:0.5706 C:0.6863 C:0.2804 C:0.6992 C:0.7056 C:0.6777 C:0.6334
2 44152180 44152180 G/A + rs4952694 A:0.9017 A:0.6398 A:0.9663 A:0.4891 A:0.6789 A:0.8207 A:0.4759 A:0.5881 A:0.8436 A:0.6603 A:0.5199 A:0.9648 A:0.4714 A:0.4931 A:0.5804 A:0.6215
2 73302892 73302892 -/GCCT + rs3832033 GCCT:0.5303 GCCT:0.2767 GCCT:0.5357 GCCT:0.2177 GCCT:0.2812 GCCT:0.4488 GCCT:0.1983 GCCT:0.219 GCCT:0.4651 GCCT:0.1816 GCCT:0.1818 GCCT:0.4731 GCCT:0.223 GCCT:0.1538 GCCT:0.2093 GCCT:0.2284
2 90121675 90121675 T/C + rs7371197 C:0.5893 C:0.9121 C:0.9236 C:0.9056 C:0.9438 C:0.5478 C:0.873 C:0.906 C:0.6271 C:0.9504 C:0.8951 C:0.9384 C:0.9398 C:0.9118 C:0.9055 C:0.929
2 98736225 98736225 C/T + rs2305355 T:0.2496 T:0.0548 T:0.1329 T:0.0656 T:0.0787 T:0.1949 T:0.07288 T:0.07964 T:0.2135 T:0.03873 T:0.1043 T:0.1412 T:0.04164 T:0.07056 T:0.07195 T:0.07692
2 98916568 98916568 T/G + rs6731704 G:0 G:0 G:0 G:0.003 G:0 G:0.0002655 G:0.002913 G:0.001389 G:0.0003271 G:0.0006254 G:0 G:0 G:0.0006726 G:0.002606 G:0.0009124 G:0.0001625
2 102626187 102626187 G/A + rs2230401 A:0.1551 A:0.013 A:0 A:0.001 A:0 A:0.1151 A:0.001395 A:0.009919 A:0.1242 A:0.009827 A:0.001624 A:0 A:0 A:0.001119 A:0.01058 A:0.0003898
20 9510263 9510263 T/A + rs2232268 A:0.7035 A:0.2089 A:0.2123 A:0.1879 A:0.2076 A:0.6158 A:0.1928 A:0.2117 A:0.6351 A:0.1222 A:0.2371 A:0.2272 A:0.1707 A:0.1847 A:0.1935 A:0.2104
20 10026357 10026357 T/C + rs685723 C:0.3048 C:0.1628 C:0.006 C:0.2763 C:0.2495 C:0.296 C:0.2806 C:0.2281 C:0.2972 C:0.1414 C:0.2604 C:0.008469 C:0.2064 C:0.2723 C:0.242 C:0.254
20 10036202 10036202 G/A + rs6087119 A:0.0461 A:0.0548 A:0 A:0.1113 A:0.0358 A:0.06128 A:0.1264 A:0.08959 A:0.06242 A:0.04623 A:0.08809 A:0.0002321 A:0.1391 A:0.1224 A:0.09532 A:0.04529
20 25276297 25276297 G/A + rs1130694 A:0.3729 A:0.3314 A:0.9107 A:0.4264 A:0.5031 A:0.3636 A:0.4256 A:0.4439 A:0.3641 A:0.3044 A:0.3111 A:0.9244 A:0.449 A:0.4289 A:0.4207 A:0.464
20 30753220 30753220 G/A + rs7273342 A:0.0378 A:0.0058 A:0 A:0 A:0 A:0.03041 A:0.0001163 A:0.00212 A:0.03045 A:0.00134 A:0.0002031 A:0 A:0 A:4.476e-05 A:0.0003646 A:6.497e-05
20 31022469 31022469 G/A + rs3746609 A:0.003 A:0.0231 A:0.0585 A:0.005 A:0.002 A:0.001953 A:0.002077 A:0.01924 A:0.001291 A:0.0815 A:0.00738 A:0.06438 A:0.009578 A:0.00248 A:0.01802 A:0.002401
20 36767986 36767986 C/T + rs2229470 T:0.0008 T:0.0086 T:0 T:0.007 T:0.0072 T:0.001362 T:0.009535 T:0.008691 T:0.00183 T:0.005568 T:0.03888 T:0.000116 T:0.002018 T:0.009375 T:0.01203 T:0.01241
20 37654024 37654024 T/G + rs3752301 G:0.09 G:0.1513 G:0.1468 G:0.2177 G:0.2342 G:0.116 G:0.1971 G:0.1768 G:0.1106 G:0.1372 G:0.2941 G:0.1683 G:0.1308 G:0.1924 G:0.1947 G:0.1936
20 40049099 40049099 G/A + rs6029672 A:0.475 A:0.2363 A:0.3631 A:0.2565 A:0.2352 A:0.4433 A:0.3085 A:0.3079 A:0.4526 A:0.2911 A:0.3116 A:0.4027 A:0.209 A:0.2973 A:0.2939 A:0.2545
20 44403108 44403108 C/G + rs6032525 G:0.0605 G:0.0086 G:0 G:0 G:0.002 G:0.05878 G:0.0001163 G:0.004232 G:0.06085 G:0.002502 G:0 G:0 G:0 G:8.989e-05 G:0.001826 G:0.0001624
21 15748157 15748157 G/A + rs9305297 A:0.0711 A:0.0043 A:0 A:0.001 A:0 A:0.04176 A:0.0002326 A:0.003516 A:0.05097 A:0.001676 A:0 A:0 A:0 A:0.0001088 A:0.001855 A:0
21 34837723 34837723 C/T + rs8130421 T:0.0204 T:0.0029 T:0 T:0 T:0 T:0.02065 T:0 T:0.001681 T:0.02368 T:0.001091 T:0 T:5.839e-05 T:0 T:1.812e-05 T:0.001117 T:0.0001022
21 41032804 41032804 G/A + rs734413 A:0.3064 A:0.7147 A:0.6468 A:0.8151 A:0.7965 A:0.3768 A:0.809 A:0.7593 A:0.3548 A:0.7645 A:0.7834 A:0.6582 A:0.8117 A:0.8044 A:0.7806 A:0.7996
21 44163839 44163839 A/G + rs6586344 G:0.9977 G:0.9885 G:0.9494 G:0.9851 G:0.8875 G:0.9959 G:0.983 G:0.9698 G:0.9959 G:0.9908 G:0.962 G:0.9561 G:0.9922 G:0.9839 G:0.9727 G:0.8933
21 44841052 44841052 T/C + rs2006112 C:0.7481 C:0.232 C:0.2718 C:0.3062 C:0.3609 C:0.6757 C:0.294 C:0.3288 C:0.6946 C:0.2016 C:0.3725 C:0.2537 C:0.3384 C:0.3102 C:0.3313 C:0.3762
21 46047971 46047971 T/C + rs463217 C:0.8714 C:0.7493 C:0.7103 C:0.7763 C:0.7832 C:0.8515 C:0.7307 C:0.7539 C:0.8465 C:0.7788 C:0.7947 C:0.7098 C:0.7749 C:0.7318 C:0.7554 C:0.7572
21 46896294 46896294 G/A + rs2230686 A:0.23 A:0.0303 A:0.0794 A:0.007 A:0.0695 A:0.1961 A:0.007114 A:0.02946 A:0.2199 A:0.01318 A:0.006415 A:0.07195 A:0.001632 A:0.005267 A:0.01782 A:0.04664
22 18910756 18910756 C/G + rs2008878 G:0.947 G:0.6787 G:0.9702 G:0.5567 G:0.7117 G:0.9104 G:0.5732 G:0.6391 G:0.9029 G:0.6248 G:0.6735 G:0.9665 G:0.4886 G:0.5716 G:0.6198 G:0.6818
22 20132972 20132972 A/G + rs175179 G:0.3986 G:0.5317 G:0.6875 G:0.3807 G:0.5194 G:0.3882 G:0.363 G:0.4561 G:0.4145 G:0.615 G:0.4385 G:0.6997 G:0.4487 G:0.3987 G:0.4599 G:0.5366
22 21344002 21344002 A/G + rs7410444 G:0.8676 G:0.9207 G:0.9861 G:0.8499 G:0.7965 G:0.8548 G:0.8373 G:0.004019 G:0.004028 G:0.0004258 G:0.0004154 G:0.001388 G:0.004883 G:0.001513 G:0.003399 G:0.01908
22 29656389 29656389 G/T + rs2231398 T:0.1006 T:0.245 T:0.0397 T:0.3608 T:0.1554 T:0.1432 T:0.3775 T:0.2758 T:0.1254 T:0.1894 T:0.4651 T:0.03356 T:0.3001 T:0.3675 T:0.3354 T:0.1643
22 38934515 38934515 C/G + rs4987164 G:0.0991 G:0.0159 G:0.1171 G:0.0159 G:0.0378 G:0.09623 G:0.01977 G:0.0297 G:0.09458 G:0.01171 G:0.02543 G:0.1043 G:0.006144 G:0.02084 G:0.02155 G:0.02736
22 39777680 39777680 C/T + rs8135738 T:0.0008 T:0.0144 T:0 T:0.0099 T:0.0174 T:0.002726 T:0.01769 T:0.01276 T:0.003059 T:0.008342 T:0.03073 T:0 T:0.004228 T:0.01541 T:0.01229 T:0.02067
22 41739370 41739370 T/C + rs1883828 C:0.7685 C:0.6988 C:0.7351 C:0.5159 C:0.636 C:0.7388 C:0.502 C:0.5942 C:0.7476 C:0.7781 C:0.5485 C:0.7385 C:0.5647 C:0.4978 C:0.5961 C:0.6113
22 42463814 42463814 C/T + rs133369 T:0.5726 T:0.6772 T:0.8532 T:0.668 T:0.6452 T:0.5976 T:0.6656 T:0.6747 T:0.5871 T:0.7415 T:0.6031 T:0.8553 T:0.6128 T:0.6719 T:0.6551 T:0.6231
22 42908870 42908870 T/G + rs9745 G:0.0582 G:0.1787 G:0.003 G:0.1461 G:0.0327 G:0.07876 G:0.1349 G:0.1321 G:0.07341 G:0.2526 G:0.05793 G:0.0008353 G:0.2035 G:0.1401 G:0.1399 G:0.04716
22 43459855 43459855 G/A + rs2272869 A:0.0008 A:0.0836 A:0.1508 A:0 A:0.0072 A:0.001589 A:0.0005814 A:0.02548 A:0.001176 A:0.1125 A:0 A:0.1335 A:4.485e-05 A:0.0001075 A:0.01751 A:0.002144
3 4562667 4562667 T/G + rs1038639 G:0.0628 G:0.4986 G:0.7044 G:0.5716 G:0.7638 G:0.1387 G:0.574 G:0.5583 G:0.1254 G:0.4754 G:0.6497 G:0.721 G:0.5028 G:0.573 G:0.5808 G:0.726
3 9885615 9885615 G/C + rs2290306 C:0.115 C:0.1037 C:0.0833 C:0.0567 C:0.1074 C:0.1011 C:0.06409 C:0.08441 C:0.1097 C:0.1321 C:0.08277 C:0.08804 C:0.04974 C:0.06877 C:0.08391 C:0.09503
3 11301791 11301791 A/G + rs2067468 G:0.0219 G:0.0115 G:0 G:0.0328 G:0.044 G:0.02202 G:0.03942 G:0.03058 G:0.02341 G:0.02148 G:0.02358 G:0.0002319 G:0.02736 G:0.03583 G:0.03522 G:0.04581
3 14157890 14157890 C/G + rs2607747 G:0.4402 G:0.1311 G:0 G:0.1322 G:0.089 G:0.3577 G:0.1378 G:0.1228 G:0.3732 G:0.08583 G:0.09777 G:0.0003479 G:0.08033 G:0.1383 G:0.1286 G:0.08902
3 19556927 19556927 T/A + rs4130926 A:0.0696 A:0.062 A:0.0714 A:0.0487 A:0.089 A:0.06971 A:0.04803 A:0.05269 A:0.06915 A:0.0448 A:0.04569 A:0.06944 A:0.05033 A:0.04482 A:0.04693 A:0.08151
3 38049997 38049997 G/A + rs2230477 A:0.0832 A:0.0029 A:0 A:0 A:0.001 A:0.064 A:0.0004651 A:0.004534 A:0.06282 A:0.002591 A:0.002336 A:0 A:0 A:0.000206 A:0.003101 A:0.0001624
3 39543794 39543794 C/T + rs2292504 T:0.0144 T:0.0029 T:0.0337 T:0.001 T:0.0061 T:0.01476 T:0.0003489 T:0.005285 T:0.01633 T:0.001642 T:0 T:0.0318 T:0.000181 T:0.0004986 T:0.003133 T:0.01195
3 42738297 42738297 G/A + rs2290973 A:0.2027 A:0.5965 A:0.6627 A:0.5378 A:0.5112 A:0.2742 A:0.5356 A:0.5317 A:0.2511 A:0.5976 A:0.478 A:0.6802 A:0.5694 A:0.5358 A:0.5405 A:0.4865
3 45137087 45137087 T/G + rs2037358 G:0.1293 G:0.6484 G:0.7966 G:0.4861 G:0.6452 G:0.1877 G:0.4817 G:0.5388 G:0.1739 G:0.7132 G:0.5083 G:0.8151 G:0.4779 G:0.4776 G:0.5253 G:0.6476
3 46488854 46488854 T/C + rs2239692 C:0.1982 C:0.1138 C:0.1984 C:0.1183 C:0.2658 C:0.1879 C:0.1139 C:0.1457 C:0.1834 C:0.09869 C:0.1276 C:0.2177 C:0.1771 C:0.1098 C:0.1207 C:0.2557
4 1360238 1360238 G/- + rs3216721 -:0.0779 -:0.0418 -:0.1885 -:0.0169 -:0.0992 -:0.07571 -:0.01527 -:0.04599 -:0.07879 -:0.02657 -:0.03623 -:0.2054 -:0.0172 -:0.01536 -:0.04301 -:0.09756
4 3487435 3487435 G/A + rs2006802 A:0.2572 A:0.036 A:0.0069 A:0.008 A:0.0051 A:0.1834 A:0.004496 A:0.01664 A:0.2069 A:0.01355 A:0.0127 A:0.004996 A:6.817e-05 A:0.003156 A:0.01593 A:0.00621
4 30725696 30725696 G/A + rs1047012 A:0.0219 A:0.1398 A:0.3105 A:0.1133 A:0.2474 A:0.03586 A:0.1227 A:0.1535 A:0.03177 A:0.1938 A:0.1074 A:0.2938 A:0.1366 A:0.1254 A:0.143 A:0.2226
4 37904082 37904082 T/C + rs2279026 C:0.6377 C:0.3444 C:0.373 C:0.2038 C:0.3701 C:0.5647 C:0.216 C:0.285 C:0.5848 C:0.3408 C:0.2028 C:0.365 C:0.2628 C:0.2164 C:0.2665 C:0.3248
4 40428010 40428010 T/C + rs278981 C:0.6725 C:0.7767 C:0.5784 C:0.7843 C:0.682 C:0.7025 C:0.7643 C:0.7494 C:0.6931 C:0.8158 C:0.8158 C:0.5877 C:0.7675 C:0.7681 C:0.7646 C:0.691
4 57301681 57301681 G/A + rs6843073 A:0.0817 A:0.111 A:0.0724 A:0.0328 A:0.1186 A:0.07425 A:0.03547 A:0.06059 A:0.07371 A:0.143 A:0.02767 A:0.05215 A:0.02511 A:0.03702 A:0.05395 A:0.09114
4 57349282 57349282 A/G + rs7700034 G:0.5885 G:0.5274 G:0.6667 G:0.329 G:0.4877 G:0.5309 G:0.2911 G:0.3991 G:0.541 G:0.5919 G:0.2144 G:0.666 G:0.4047 G:0.293 G:0.3792 G:0.4189
4 57384939 57384939 T/C + rs7672037 C:0.0045 C:0 C:0 C:0 C:0 C:0.002106 C:0 C:0.0002153 C:0.003336 C:2.978e-05 C:0 C:0 C:0 C:0 C:0.0001825 C:0
4 71892513 71892513 A/T + rs1486271 T:0.3328 T:0.9049 T:0.9603 T:0.9612 T:0.9591 T:0.4508 T:0.9665 T:0.9257 T:0.3933 T:0.9442 T:0.9633 T:0.9488 T:0.972 T:0.968 T:0.9457 T:0.9623
4 83857108 83857108 A/G + rs6535413 G:0.8139 G:0.3487 G:0.7986 G:0.4304 G:0.6493 G:0.09169 G:0.0314 G:0.4776 G:0.7669 G:0.2874 G:0.5243 G:0.8005 G:0.3568 G:0.4277 G:0.4535 G:0.6164
5 5200281 5200281 C/T + rs6555335 T:0.5651 T:0.8314 T:0.8026 T:0.7396 T:0.772 T:0.6145 T:0.7477 T:0.7577 T:0.5859 T:0.8633 T:0.6692 T:0.8067 T:0.7446 T:0.7478 T:0.7463 T:0.7764
5 54572163 54572163 T/C + rs3761764 C:0.2943 C:0.0965 C:0.2877 C:0.0994 C:0.0798 C:0.2234 C:0.1012 C:0.1103 C:0.236 C:0.05552 C:0.09445 C:0.2665 C:0.1002 C:0.09746 C:0.1104 C:0.07741
5 54646894 54646894 A/G + rs9885400 G:0.292 G:0.098 G:0.2867 G:0.0994 G:0.0787 G:0.2189 G:0.1003 G:0.1202 G:0.249 G:0.06674 G:0.1013 G:0.2784 G:0.1032 G:0.102 G:0.1228 G:0.08718
5 66479896 66479896 C/T + rs5744525 T:0.003 T:0.0014 T:0 T:0 T:0 T:0.00749 T:0.0001163 T:0.0003101 T:0.004515 T:0.0001192 T:0 T:0 T:0 T:1.808e-05 T:0.0001827 T:0
5 72347223 72347223 A/G + rs7712838 G:0.2897 G:0.1859 G:0.4782 G:0.3181 G:0.2035 G:0.266 G:0.2721 G:0.2724 G:0.2812 G:0.1701 G:0.2512 G:0.4475 G:0.2855 G:0.3055 G:0.2871 G:0.2082
5 76129241 76129241 G/A + rs2243062 A:0.0431 A:0.0043 A:0 A:0 A:0 A:0.04108 A:0.0004651 A:0.003338 A:0.0462 A:0.002353 A:0 A:0 A:0 A:0.0001612 A:0.002553 A:0.0001299
5 76249536 76249536 C/G + rs1715771 G:0.0832 G:0.5115 G:0.4821 G:0.4662 G:0.6759 G:0.1413 G:0.457 G:0.4988 G:0.1407 G:0.5883 G:0.3543 G:0.4811 G:0.479 G:0.5005 G:0.4458 G:0.6244
5 79734297 79734297 T/C + rs259028 C:0.6944 C:0.9294 C:0.9196 C:0.9334 C:0.865 C:0.7199 C:0.9434 C:0.9146 C:0.7066 C:0.9537 C:0.9463 C:0.9311 C:0.9222 C:0.9419 C:0.9198 C:0.851
5 118691832 118691832 C/T + rs3797342 T:0.2156 T:0.2291 T:0.0516 T:0.2048 T:0.1585 T:0.2092 T:0.2058 T:0.1933 T:0.2067 T:0.1981 T:0.2199 T:0.05316 T:0.2618 T:0.2033 T:0.2091 T:0.1599
5 126746277 126746277 C/T + rs31483 T:0.9977 T:0.9568 T:1 T:0.9324 T:0.9673 T:0.9821 T:0.9199 T:0.9443 T:0.9872 T:0.9665 T:0.9508 T:0.9996 T:0.9251 T:0.9202 T:0.9377 T:0.9645
6 6174786 6174786 G/C + rs380058 C:0.9992 C:0.8732 C:0.9454 C:0.999 C:0.999 C:0.9986 C:0.999 C:0.9698 C:0.998 C:0.8243 C:0.9998 C:0.9449 C:0.9888 C:0.9992 C:0.9768 C:0.9988
6 10621547 10621547 C/T + rs9460944 T:0.3094 T:0.1239 T:0.002 T:0.0378 T:0.0348 T:0.2773 T:0.04465 T:0.06628 T:0.2793 T:0.1013 T:0.02144 T:0.006647 T:0.05492 T:0.04957 T:0.06473 T:0.04419
6 10751379 10751379 C/T + rs1046427 T:0.0779 T:0.0115 T:0 T:0.004 T:0.002 T:0.05949 T:0.003256 T:0.007323 T:0.06489 T:0.006522 T:0.0132 T:0 T:0.0002691 T:0.003367 T:0.00857 T:0.00104
6 13288669 13288669 C/T + rs202037 T:0.034 T:0.2666 T:0.0179 T:0.0706 T:0.1043 T:0.04852 T:0.07308 T:0.1181 T:0.0376 T:0.3791 T:0.07682 T:0.02248 T:0.06732 T:0.07741 T:0.1143 T:0.09575
6 17764731 17764731 G/C + rs6459569 C:0.3979 C:0.098 C:0.0635 C:0.0696 C:0.1043 C:0.3399 C:0.08335 C:0.09957 C:0.3698 C:0.08705 C:0.1001 C:0.06321 C:0.04991 C:0.08051 C:0.09036 C:0.106
6 25420344 25420344 C/G + rs913455 G:0.8631 G:0.9625 G:1 G:0.9443 G:1 G:0.8828 G:0.9343 G:0.9545 G:0.873 G:0.9738 G:0.9752 G:1 G:0.9688 G:0.9357 G:0.9589 G:0.9989
6 28359073 28359073 A/G + rs2232432 G:0.0197 G:0.0086 G:0 G:0.008 G:0.001 G:0.01084 G:0.007542 G:0.00779 G:0.01578 G:0.0103 G:0.0451 G:9.731e-05 G:0.0004772 G:0.005822 G:0.01498 G:0.001533
6 30596135 30596135 A/G + rs6904236 G:0.1157 G:0.0346 G:0.0258 G:0.0199 G:0.0082 G:0.1019 G:0.02049 G:0.0239 G:0.1119 G:0.02907 G:0.04589 G:0.01927 G:0.001391 G:0.01804 G:0.02847 G:0.009295
6 30893728 30893728 C/T + rs1043483 T:0.4622 T:0.696 T:0.8671 T:0.7137 T:0.8783 T:0.5096 T:0.7329 T:0.7541 T:0.4904 T:0.7472 T:0.7732 T:0.8788 T:0.7715 T:0.737 T:0.7282 T:0.8607
6 31116502 31116502 C/T + rs2027937 T:0.0378 T:0.1023 T:0.0417 T:0.0795 T:0.1237 T:0.05199 T:0.07254 T:0.08356 T:0.04293 T:0.09307 T:0.04799 T:0.04116 T:0.09562 T:0.07643 T:0.0936 T:0.1417
7 1533556 1533556 T/C + rs2251231 C:0.5908 C:0.3876 C:0.7123 C:0.4493 C:0.499 C:0.5076 C:0.3782 C:0.4539 C:0.563 C:0.3725 C:0.3928 C:0.6797 C:0.4734 C:0.4134 C:0.4273 C:0.5171
7 22233674 22233674 T/C + rs6965394 C:0.0318 C:0.0432 C:0.0089 C:0.1044 C:0.0726 C:0.03071 C:0.08132 C:0.07451 C:0.03513 C:0.04589 C:0.1164 C:0.006232 C:0.1215 C:0.08549 C:0.08833 C:0.07246
7 37901581 37901581 G/A + rs2249451 A:0.916 A:0.8804 A:0.6508 A:0.9453 A:0.7791 A:0.9165 A:0.9441 A:0.8854 A:0.9186 A:0.8117 A:0.9631 A:0.6892 A:0.8614 A:0.9464 A:0.9064 A:0.8275
7 44259706 44259706 G/A + rs1065359 A:0.1596 A:0.3213 A:0.3879 A:0.4085 A:0.4724 A:0.2207 A:0.4174 A:0.4064 A:0.2029 A:0.3733 A:0.4501 A:0.3447 A:0.508 A:0.4229 A:0.3996 A:0.4369
7 75145495 75145495 G/A + rs794377 A:0.0189 A:0.598 A:0.6647 A:0.4443 A:0.498 A:0.07944 A:0.3692 A:0.4534 A:0.07387 A:0.6468 A:0.3847 A:0.658 A:0.5248 A:0.4095 A:0.4586 A:0.4381
7 76891485 76891485 A/G + rs1109968 G:0.1082 G:0.2522 G:0.4206 G:0.0467 G:0.1237 G:0.0901 G:0.04651 G:0.1257 G:0.09338 G:0.3139 G:0.06979 G:0.4257 G:0.06343 G:0.0472 G:0.1082 G:0.1203
7 92147374 92147374 C/T + rs2066743 T:0.0113 T:0.098 T:0.2758 T:0.0189 T:0.0317 T:0.01657 T:0.02256 T:0.05806 T:0.0184 T:0.1557 T:0.06927 T:0.2634 T:0.004576 T:0.01968 T:0.04403 T:0.03311
7 99709314 99709314 C/T + rs4134906 T:0.0741 T:0.0043 T:0 T:0 T:0 T:0.06083 T:0.0001163 T:0.004779 T:0.06454 T:0.003068 T:0 T:0 T:0 T:0.0004929 T:0.003832 T:0.0002924
7 100678740 100678740 T/C + rs4269454 C:0.351 C:0.1023 C:0.252 C:0.1213 C:0.2883 C:0.3059 C:0.1162 C:0.1585 C:0.318 C:0.1334 C:0.1333 C:0.2184 C:0.1198 C:0.1143 C:0.1449 C:0.2713
7 103824667 103824667 A/C + rs7807228 C:0.0378 C:0.0029 C:0 C:0.001 C:0 C:0.02387 C:0 C:0.002759 C:0.0347 C:0.002317 C:0 C:0 C:0 C:0.0001282 C:0.002841 C:0
8 2836102 2836102 C/G + rs615578 G:0.77 G:0.6628 G:0.8056 G:0.8231 G:0.7301 G:0.7996 G:0.817 G:0.7929 G:0.7901 G:0.6803 G:0.8093 G:0.8178 G:0.871 G:0.8195 G:0.7841 G:0.7528
8 11643459 11643459 C/T + rs8191663 T:0.3631 T:0.1859 T:0.2738 T:0.2634 T:0.1575 T:0.3196 T:0.2334 T:0.2272 T:0.3236 T:0.1579 T:0.2918 T:0.2832 T:0.2255 T:0.2312 T:0.2461 T:0.1871
8 11710974 11710974 G/A + rs1137063 A:0.0083 A:0.0562 A:0 A:0.0905 A:0.0174 A:0.01566 A:0.08477 A:0.06199 A:0.01189 A:0.04385 A:0.07011 A:0.0004076 A:0.08318 A:0.08955 A:0.08121 A:0.02134
8 25257531 25257531 A/G + rs3763520 G:0.2042 G:0.6225 G:0.6339 G:0.6332 G:0.5849 G:0.2812 G:0.632 G:0.6005 G:0.2618 G:0.6259 G:0.6099 G:0.6159 G:0.6389 G:0.6268 G:0.6025 G:0.6047
8 27348788 27348788 C/A + rs2234912 A:0.0008 A:0.0058 A:0 A:0.003 A:0.001 A:0.001664 A:0.004952 A:0.004433 A:0.001356 A:0.002689 A:0 A:0 A:0.009664 A:0.006982 A:0.005909 A:0.0009527
8 37728019 37728019 T/G + rs7820872 G:0.8071 G:0.7738 G:0.9752 G:0.7416 G:0.6922 G:0.7658 G:0.719 G:0.7563 G:0.7767 G:0.8317 G:0.6876 G:0.9766 G:0.7736 G:0.7173 G:0.7428 G:0.6937
8 42177163 42177163 G/A + rs2272736 A:0 A:0 A:0.1181 A:0 A:0 A:0.001362 A:0.0001163 A:0.008529 A:0.000875 A:0.0001881 A:0.0002416 A:0.1136 A:4.624e-05 A:0.0001139 A:0.003919 A:0.002977
8 72942210 72942210 G/T + rs3824151 T:0.6732 T:0.7536 T:0.6935 T:0.5726 T:0.7127 T:0.6534 T:0.5905 T:0.643 T:0.6513 T:0.8075 T:0.6566 T:0.6849 T:0.5596 T:0.5935 T:0.6474 T:0.6697
8 121228679 121228679 A/C + rs4870723 C:0.6286 C:0.5749 C:0.6687 C:0.4702 C:0.591 C:0.5756 C:0.4855 C:0.5128 C:0.5784 C:0.5955 C:0.486 C:0.6386 C:0.4033 C:0.474 C:0.5066 C:0.5492
8 124448736 124448736 A/G + rs7014678 G:0.3427 G:0.3919 G:0.5863 G:0.3111 G:0.4346 G:0.3631 G:0.3522 G:0.3813 G:0.3579 G:0.4164 G:0.3864 G:0.5546 G:0.3772 G:0.3412 G:0.3779 G:0.4051
9 6012968 6012968 A/G + rs1411949 G:0.3321 G:0.0375 G:0.0714 G:0.002 G:0.1554 G:0.2667 G:0.002326 G:0.04265 G:0.2644 G:0.01525 G:0.009748 G:0.0534 G:0.003818 G:0.002508 G:0.02647 G:0.1435
9 35957728 35957728 C/T + rs2233563 T:0.0915 T:0.0014 T:0 T:0 T:0 T:0.07966 T:0.0003488 T:0.006109 T:0.08429 T:0.004378 T:0 T:0 T:0 T:0.0003135 T:0.005108 T:0.00013
9 73442724 73442724 C/T + rs1034538 T:0.9569 T:0.8026 T:0.5427 T:0.7763 T:0.7025 T:0.9094 T:0.7787 T:0.7591 T:0.9156 T:0.789 T:0.8266 T:0.5305 T:0.7907 T:0.7675 T:0.77 T:0.7006
9 79320640 79320640 C/T + rs512110 T:0.739 T:0.8228 T:0.9117 T:0.7734 T:0.7628 T:0.7497 T:0.7627 T:0.8046 T:0.7508 T:0.8735 T:0.8632 T:0.9381 T:0.8201 T:0.7727 T:0.8108 T:0.7655
9 93641151 93641151 T/C + rs2306042 C:0.0008 C:0.0072 C:0.2321 C:0.0099 C:0.1135 C:0.003631 C:0.01337 C:0.03756 C:0.003332 C:0.00804 C:0.01361 C:0.2343 C:0.007445 C:0.01591 C:0.02243 C:0.08729
9 107360851 107360851 T/C + rs1523678 C:0.6475 C:0.2277 C:0.5754 C:0.1789 C:0.4192 C:0.5849 C:0.1941 C:0.2773 C:0.598 C:0.2734 C:0.1179 C:0.5749 C:0.2451 C:0.196 C:0.2292 C:0.3325
9 115598447 115598447 C/T + rs6477961 T:0.5961 T:0.7896 T:0.7877 T:0.7922 T:0.8609 T:0.6372 T:0.8045 T:0.7917 T:0.6236 T:0.7355 T:0.836 T:0.8211 T:0.8378 T:0.7988 T:0.8024 T:0.8584
9 115920099 115920099 T/C + rs4344146 C:0.671 C:0.6095 C:0.5952 C:0.5586 C:0.5511 C:0.6823 C:0.5765 C:0.5593 C:0.6732 C:0.5201 C:0.5958 C:0.6084 C:0.4493 C:0.5647 C:0.5477 C:0.5699
9 116857656 116857656 T/C + rs7036734 C:0.5635 C:0.5576 C:0.7083 C:0.5 C:0.7014 C:0.5423 C:0.5034 C:0.5687 C:0.5564 C:0.603 C:0.5812 C:0.6875 C:0.5282 C:0.5089 C:0.5552 C:0.6821
9 116931666 116931666 A/T + rs2567705 T:0.3654 T:0.2795 T:0.1151 T:0.4215 T:0.2454 T:0.3777 T:0.3916 T:0.3315 T:0.3777 T:0.2116 T:0.5965 T:0.1303 T:0.3147 T:0.3835 T:0.3709 T:0.2843
X 2700202 2700202 A/G + rs5939320 G:0.3988 G:0.2882 G:0.5092 G:0.3042 G:0.2242 G:0.4141 G:0.3091 G:0.3055 G:0.4071 G:0.2307 G:0.2192 G:0.5224 G:0.2655 G:0.3108 G:0.3089 G:0.2433
X 35990125 35990125 A/- + rs5902124 -:0.7348 -:0.3531 -:0.3037 -:0.0822 -:0.2744 -:0.6843 -:0.06662 -:0.1842 -:0.7141 -:0.2644 -:0.1024 -:0.3099 -:0.116 -:0.06545 -:0.1741 -:0.246
X 38080703 38080703 G/T + rs7877137 T:0.0289 T:0.2023 T:0.055 T:0.4282 T:0.2577 T:0.0811 T:0.4062 T:0.2922 T:0.05953 T:0.1429 T:0.4813 T:0.04497 T:0.3697 T:0.4016 T:0.3561 T:0.2669
X 53457622 53457622 T/C + rs1264013 C:0.1675 C:0.521 C:0.538 C:0.6031 C:0.5 C:0.2068 C:0.5871 C:0.5379 C:0.2003 C:0.5287 C:0.4898 C:0.5304 C:0.6026 C:0.595 C:0.5394 C:0.4994
X 54117740 54117740 C/A + rs2495797 A:0.99 A:0.9981 A:1 A:1 A:1 A:0.988 A:0.9999 A:0.9987 A:0.9838 A:0.999 A:1 A:1 A:1 A:0.9999 A:1 A:1
X 71521798 71521798 T/C + rs1804686 C:0.341 C:0.0286 C:0 C:0.0026 C:0 C:0.2615 C:0.001189 C:0.02136 C:0.2707 C:0.01419 C:0 C:0 C:0 C:0.0007927 C:0.01291 C:0.001681
X 85134098 85134098 G/C + rs5968705 C:0.3121 C:0.0439 C:0.0825 C:0 C:0.0223 C:0.2579 C:0.002245 C:0.03069 C:0.2758 C:0.01551 C:0.006841 C:0.08014 C:0.0002295 C:0.001829 C:0.02374 C:0.009118
X 88009322 88009322 A/G + rs5984613 G:0.9741 G:0.9141 G:0.9228 G:0.859 G:0.8774 G:0.9679 G:0.8728 G:0.8868 G:0.9704 G:0.9033 G:0.8922 G:0.9208 G:0.8632 G:0.8718 G:0.8777 G:0.8524
X 100400023 100400023 A/G + rs2295376 G:0.9611 G:0.5992 G:0.3966 G:0.7637 G:0.6741 G:0.9233 G:0.7491 G:0.6946 G:0.9381 G:0.5004 G:0.7837 G:0.4252 G:0.6976 G:0.7452 G:0.7195 G:0.6982
X 111698613 111698613 T/C + rs7053563 C:0.5005 C:0.292 C:0.3455 C:0.2089 C:0.273 C:0.4326 C:0.2012 C:0.2621 C:0.4421 C:0.302 C:0.185 C:0.3702 C:0.2804 C:0.2042 C:0.2493 C:0.2758
Y 14832620 14832620 G/T + rs7067496 T:0.931 T:0.1824 T:0.1311 T:0.025 T:0.0692 T:0.69 T:0.05716 T:0.08915 T:0.7451 T:0.09489 T:0.2131 T:0.08762 T:0.002907 T:0.05014 T:0.105 T:0.0332
Y 14851554 14851554 T/C + rs2032599 C:0.0188 C:0.0059 C:0 C:0 C:0 C:0.02452 C:0 C:0.001198 C:0.02299 C:0.0007298 C:0 C:0 C:0 C:6.402e-05 C:0.00142 C:8.686e-05
Y 14869076 14869076 C/T + rs2032601 T:0.0094 T:0.0059 T:0 T:0 T:0 T:0.01576 T:0 T:0.0008214 T:0.01404 T:0.0007344 T:0 T:0.0002595 T:0 T:0.0001098 T:0.002271 T:9.159e-05
Y 14888783 14888783 A/C + rs2032600 C:0 C:0.0412 C:0 C:0.0042 C:0 C:0.001751 C:0.004274 C:0.004949 C:0.00131 C:0.02136 C:0.004317 C:0 C:0.0001729 C:0.003886 C:0.03123 C:8.61e-05
Y 14902414 14902414 G/A + rs20321 A:0.0031 A:0 A:0 A:0 A:0 A:0.003503 A:0.03579 A:0.009614 A:0.006544 A:0.0008574 A:0.000392 A:0 A:0.0001708 A:0.01933 A:0.007725 A:8.598e-05
Y 15469724 15469724 G/A + rs9341278 A:0 A:0 A:0.0328 A:0.0958 A:0.0038 A:0 A:0.005876 A:0.07612 A:0.0003937 A:0.001661 A:0.0004146 A:0.0409 A:0.6944 A:0.01712 A:0.07514 A:0.00051
Y 15591537 15591537 G/C + rs2032653 C:0.8966 C:0.1765 C:0.082 C:0.025 C:0 C:0.6567 C:0.05662 C:0.07925 C:0.7276 C:0.09543 C:0.2131 C:0.01325 C:0.002578 C:0.04961 C:0.1046 C:0.01038
Y 21869264 21869264 G/A + rs2032620 A:0.0094 A:0 A:0 A:0 A:0 A:0.02977 A:0 A:0.001341 A:0.02715 A:0.0007612 A:0 A:0.0004965 A:0 A:3.679e-05 A:0.001555 A:0
Y 21893881 21893881 A/C + rs2032672 C:0 C:0.0235 C:0 C:0.0167 C:0 C:0.001751 C:0.006944 C:0.008169 C:0.0003272 C:0.01461 C:0.03532 C:0 C:0.0001713 C:0.008068 C:0.01758 C:0.006622
};
  }
  else {
    die "Unrecognised assembly $assembly\n";
  }

  close OUT;

  return $file;
}

sub DESTROY {
  unlink($_[0]->{_lock_file}) if $_[0]->{_lock_file};
}

1;
