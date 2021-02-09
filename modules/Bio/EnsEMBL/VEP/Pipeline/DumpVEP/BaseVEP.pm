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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::BaseVEP;

use strict;
use warnings;

use File::Path qw(make_path rmtree);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub data_dir {
  my $self = shift;
  return $self->param('pipeline_dump_dir').$self->param('dir_suffix');
}

sub dump_dir {
  my $self = shift;
  return $self->param('pipeline_dump_dir').'/dumps'.$self->param('dir_suffix');
}

sub link_dir_contents {
  my ($self, $source_dir, $target_dir, $core, $var, $reg) = @_;

  opendir SOURCE, $source_dir or die $!;
  my @chrs = grep {$_ !~ /^\./ && -d $source_dir.'/'.$_} readdir(SOURCE);

  foreach my $chr(@chrs) {
    make_path($target_dir.'/'.$chr);

    opendir CHR, $source_dir.'/'.$chr;
    my @chr_contents = readdir(CHR);
    closedir CHR;

    # core files
    if($core) {
      $self->link_file($source_dir.'/'.$chr.'/'.$_, $target_dir.'/'.$chr.'/'.$_) for grep {/\d+\-\d+\.gz$/} @chr_contents;
    }

    # var files
    if($var) {
      if($var == 1) {
        $self->link_file($source_dir.'/'.$chr.'/'.$_, $target_dir.'/'.$chr.'/'.$_) for grep {/\d+\-\d+\_var.gz$/} @chr_contents;
      }
      elsif($var == 2) {
        $self->link_file($source_dir.'/'.$chr.'/all_vars.gz', $target_dir.'/'.$chr.'/all_vars.gz');
        $self->link_file($source_dir.'/'.$chr.'/all_vars.gz.csi', $target_dir.'/'.$chr.'/all_vars.gz.csi');
      }
    }

    # reg files
    if($reg) {
      $self->link_file($source_dir.'/'.$chr.'/'.$_, $target_dir.'/'.$chr.'/'.$_) for grep {/\d+\-\d+\_reg.gz$/} @chr_contents; 
    }
  }

  closedir SOURCE;
}

sub link_file {
  my ($self, $source, $target) = @_;

  $self->run_system_command(
    sprintf(
      'ln -f %s %s',
      $source,
      $target
    )
  );
}

sub get_tar_file_name {
  my ($self, $dir, $species, $assembly, $mod) = @_;
  
  return sprintf(
    '%s/%s_vep_%i_%s%s.tar.gz',
    $dir,
    $species,
    $self->param('eg_version') || $self->param('ensembl_release'),
    $assembly,
    $mod || '',
  );
}

sub species_suffix {
  my $self = shift;

  return sprintf(
    '%s/%s_%s',
    $self->param('species'),
    $self->param('eg_version') || $self->param('ensembl_release'),
    $self->param('assembly')
  );
}

sub refseq_species_suffix {
  my $self = shift;

  return sprintf(
    '%s_refseq/%s_%s',
    $self->param('species'),
    $self->param('eg_version') || $self->param('ensembl_release'),
    $self->param('assembly')
  ); 
}

sub merged_species_suffix {
  my $self = shift;

  return sprintf(
    '%s_merged/%s_%s',
    $self->param('species'),
    $self->param('eg_version') || $self->param('ensembl_release'),
    $self->param('assembly')
  ); 
}

sub run_system_command {
  my $self = shift;
  my $cmd = shift;

  my ($exit_code, $stderr, $flat_cmd) = $self->SUPER::run_system_command($cmd);
  die("Failed to run $flat_cmd: $stderr\n") unless $exit_code == 0;

  return $exit_code;
}

1;
