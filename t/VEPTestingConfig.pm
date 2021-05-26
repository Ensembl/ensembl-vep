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

package VEPTestingConfig;

use FindBin qw($Bin);
use Compress::Zlib;

our %DEFAULTS = (
  cache_root_dir => $Bin.'/testdata/cache/',
  cache_species  => 'homo_sapiens',
  cache_version  => 84,
  cache_assembly => 'GRCh38',
  cache_dir      => $Bin.'/testdata/cache/homo_sapiens/84_GRCh38',

  cache_chr      => 21,
  cache_region   => '25000001-26000000',
  cache_s        => 25,
  var_cols       => [qw(
    variation_name failed somatic start end
    allele_string strand minor_allele minor_allele_freq
    clin_sig phenotype_or_disease pubmed
    AFR AMR EAS EUR SAS AA EA
    gnomAD gnomAD_AFR gnomAD_AMR gnomAD_ASJ gnomAD_EAS gnomAD_FIN gnomAD_NFE gnomAD_OTH gnomAD_SAS
  )],

  chr_synonyms   => $Bin.'/testdata/chr_synonyms.txt',
  haplo_freqs    => $Bin.'/testdata/haplotype_frequencies.txt',
  haplo_freqs_nh => $Bin.'/testdata/haplotype_frequencies_nh.txt',

  sereal_dir     => $Bin.'/testdata/cache/sereal/homo_sapiens/84_GRCh38',
  exac_root_dir  => $Bin.'/testdata/cache/exac/',

  fasta          => $Bin.'/testdata/cache/homo_sapiens/84_GRCh38/test.fa',
  fasta_dir      => $Bin.'/testdata/fasta',

  test_ini_file  => $Bin.'/testdata/vep.ini',
  registry_file  => $Bin.'/testdata/vep.registry',

  test_vcf       => $Bin.'/testdata/input/test.vcf',
  test_vcf2      => $Bin.'/testdata/input/test2.vcf',
  test_vcf3      => $Bin.'/testdata/input/test3.vcf',
  test_vcf4      => $Bin.'/testdata/input/test4.vcf',
  test_vcf5      => $Bin.'/testdata/input/test5.vcf',
  test_vcf6      => $Bin.'/testdata/input/test6.vcf',
  test_vcf_MT    => $Bin.'/testdata/input/test_MT.vcf',
  not_ord_vcf    => $Bin.'/testdata/input/test_not_ordered.vcf',
  no_trans_vcf   => $Bin.'/testdata/input/test_no_trans.vcf',
  vr_vcf         => $Bin.'/testdata/input/idt_test.vcf',
  windows_vcf    => $Bin.'/testdata/input/windows.vcf',
  test_gzvcf     => $Bin.'/testdata/input/test.vcf.gz',
  user_file      => $Bin.'/testdata/user_file'.$$,

  custom_bed     => $Bin.'/testdata/custom/test.bed.gz',
  custom_vcf     => $Bin.'/testdata/custom/test.vcf.gz',
  custom_gff     => $Bin.'/testdata/custom/test.gff.gz',
  custom_refseq_gff => $Bin.'/testdata/custom/refseq.gff.gz',
  custom_gtf     => $Bin.'/testdata/custom/test.gtf.gz',
  custom_gtf_mt     => $Bin.'/testdata/custom/test_MT.gtf.gz',
  custom_bigwig  => $Bin.'/testdata/custom/test.bw',

  bam_edit_gff   => $Bin.'/testdata/custom/bam_edit.gff.gz',
  bam_edit_bam   => $Bin.'/testdata/custom/bam_edit.bam',

  filter_list    => $Bin.'/testdata/filter/gene_list.txt',

  header_info    => {
    time => 'test',
    vep_version => 1,
    api_version => 1,
    db_version => 1,
    custom_info => [{
      'short_name' => 'custom_test',
      'file' => 'test.vcf.gz',
      'type' => 'overlap',
    }],
  }
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my %config = %DEFAULTS;

  # initialise self
  my $self = bless \%config, $class;

  $self->unpack_fasta;

  return $self;
}

sub unpack_fasta {
  my $self = shift;

  unless(-e $self->{cache_dir}.'/test.fa') {
    my $gz = gzopen($self->{cache_dir}.'/test.fa.gz', 'rb');

    my $buffer;
    open OUT, ">".$self->{cache_dir}.'/test.fa';
    while($gz->gzread($buffer)) {
      print OUT $buffer;
    }
    close OUT;
  }
}

# returns a hashref for general testing use
sub base_testing_cfg {
  my $self = shift;

  return {
    dir           => $self->{cache_root_dir},
    species       => $self->{cache_species},
    cache_version => $self->{cache_version},
    assembly      => $self->{cache_assembly},
    offline       => 1,
    database      => 0,
    no_stats      => 1,
  }
}

# reads MultTestDB.conf for DB params
# the same file is used by Bio::EnsEMBL::Test::MultiTestDB
sub db_cfg {
  my $self = shift;

  if(!exists($self->{db_cfg})) {

    my $cfg = {};

    if(open IN, $Bin.'/MultiTestDB.conf') {
      my @lines = <IN>;
      $cfg = eval join('', @lines);
      $cfg->{password} = $cfg->{pass};
      close IN;
    }

    $self->{db_cfg} = $cfg;
  }

  return $self->{db_cfg};
}

# creates a registry file with the params from MultTestDB.conf
# caveat is that you must pass in the dbname parameter as created by Bio::EnsEMBL::Test::MultiTestDB
# AFTER you have called Bio::EnsEMBL::Test::MultiTestDB->new()
# $vep_testing_cfg->registry_file($multi->{conf}->{core}->{dbname})
sub registry_file {
  my $self = shift;
  my $dbname = shift;

  if(!exists($self->{user_registry_file})) {

    die("ERROR: No dbname given\n") unless $dbname;

    my %db_cfg = %{$self->db_cfg()};

    die("ERROR: No db config found\n") unless scalar keys %db_cfg;

    $db_cfg{dbname} = $dbname;

    my $base_file = $self->{registry_file};

    # we need to write a new version of this file
    my $reg_file = $base_file.$$;

    open IN, $base_file or die "ERROR: Could not read from $base_file\n";
    open OUT, ">$reg_file" or die "ERROR: Could not write to $reg_file\n";

    while(<IN>) {
      my $line = $_;
      $line =~ s/\_\_$_\_\_/$db_cfg{$_}/ge for keys %db_cfg;
      print OUT $line;
    }

    close IN;
    close OUT;

    $self->{user_registry_file} = $reg_file;
  }

  return $self->{user_registry_file};
}

# creates an input file for testing
sub create_input_file {
  my $self = shift;
  my $data = shift;

  open OUT, ">".$self->{user_file} or die "ERROR: Could not write to file ".$self->{user_file}."\n";

  if($data) {
    if(ref($data) eq 'ARRAY') {
      if(ref($data->[0]) eq 'ARRAY') {
        print OUT join("\t", @$_)."\n" for @$data;
      }
      else {
        print OUT join("\t", @$data)."\n";
      }
    }
    else {
      print OUT "$data\n";
    }
  }

  close OUT;

  return $self->{user_file};
}

# remove any created files
sub DESTROY {
  my $self = shift;

  for my $file_key(qw(user_file user_registry_file)) {
    unlink($self->{$file_key}) if $self->{$file_key};
  }
}
