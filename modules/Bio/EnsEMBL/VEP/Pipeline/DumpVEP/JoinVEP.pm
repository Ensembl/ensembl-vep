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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::JoinVEP;

use strict;
use warnings;
use File::Path qw(rmtree);

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::BaseVEP);

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

  my $type = $self->param('type');

  $self->$type();

  return 1;
}

sub core {
  my $self = shift;

  my $source_dir = $self->data_dir.'/'.$self->species_suffix;
  my $target_dir = $self->dump_dir.'/'.$self->species_suffix;

  # clear any previous existing dir
  rmtree($target_dir);

  my $var = $self->param('variation') ? 1 : 0;
  my $reg = $self->param('regulation') ? 1 : 0;

  # add content
  $self->link_dir_contents($source_dir, $target_dir, 1, $var, $reg);

  # create info file
  unlink($target_dir.'/info.txt');
  $self->run_system_command(
    sprintf(
      'cat %s >> %s',
      $source_dir.'/info.txt_'.$_,
      $target_dir.'/info.txt'
    )
  ) for ('core', grep {$self->param($_)} qw(variation regulation));

  # copy synonyms file
  $self->copy_synonyms($target_dir);
  
  # Generate MD5s
  $self->generate_md5s if ($self->param('species') eq 'homo_sapiens');

  # create tar
  $self->tar(
    $self->dump_dir,
    $self->param('species'),
    $self->param('assembly'),
    $self->species_suffix
  );

  # converted?
  $self->create_converted_version(
    $source_dir,
    $target_dir,
    $source_dir,
    $self->param('species'),
    $self->species_suffix
  ) if $self->param('convert') && $var;

  rmtree($target_dir);
}

sub refseq {
  return $_[0]->_base_refseq_merged('refseq');
}

sub merged {
  return $_[0]->_base_refseq_merged('merged');
}

sub _base_refseq_merged {
  my $self = shift;
  my $type = shift;

  my $suffix_method_name = $type.'_species_suffix';

  my $source_dir = $self->data_dir.'/'.$self->$suffix_method_name;
  my $target_dir = $self->dump_dir.'/'.$self->$suffix_method_name;

  # clear any previous existing dir
  rmtree($target_dir);
  
  # copy core stuff first
  $self->link_dir_contents($source_dir, $target_dir, 1, 0, 0);

  # create info file
  unlink($target_dir.'/info.txt');
  $self->run_system_command(
    sprintf(
      'cat %s >> %s',
      $source_dir.'/info.txt_core',
      $target_dir.'/info.txt'
    )
  );

  # copy synonyms file
  $self->copy_synonyms($target_dir);

  # now copy var and reg stuff
  my $var = $self->param('variation') ? 1 : 0;
  my $reg = $self->param('regulation') ? 1 : 0;

  if($var || $reg) {
    $source_dir = $self->data_dir.'/'.$self->species_suffix;

    $self->link_dir_contents($source_dir, $target_dir, 0, $var, $reg);

    $self->run_system_command(
      sprintf(
        'cat %s >> %s',
        $source_dir.'/info.txt_'.$_,
        $target_dir.'/info.txt'
      )
    ) for grep {$self->param($_)} qw(variation regulation);
  }

  $self->tar(
    $self->dump_dir,
    $self->param('species').'_'.$type,
    $self->param('assembly'),
    $self->$suffix_method_name,    
  );

  # converted?
  $self->create_converted_version(
    $source_dir,
    $target_dir,
    $self->data_dir.'/'.$self->$suffix_method_name,
    $self->param('species').'_'.$type,
    $self->$suffix_method_name
  ) if $self->param('convert') && $var;

  rmtree($target_dir);
}

sub clear_var_files {
  my ($self, $dir) = @_;

  opendir DIR, $dir;
  my @chrs = grep {$_ !~ /^\./ && -d $dir.'/'.$_} readdir(DIR);

  foreach my $chr(@chrs) {
    $self->run_system_command(sprintf('rm -f %s/%s/*_var.gz', $dir, $chr));
  }

  closedir DIR;
}

sub create_converted_version {
  my ($self, $source_dir, $target_dir, $info_dir, $species, $suffix) = @_;

  # clear the per-MB var files
  $self->clear_var_files($target_dir);

  # now copy the all_vars.gz files
  $self->link_dir_contents($source_dir, $target_dir, 0, 2, 0);

  # create info file
  unlink($target_dir.'/info.txt');

  # $info_dir is the core-type info.txt, this might be different for _refseq and _merged caches
  $self->run_system_command(
    sprintf(
      'cat %s >> %s',
      $info_dir.'/info.txt_core',
      $target_dir.'/info.txt'
    )
  );

  # get the other info bits from the core dir
  $self->run_system_command(
    sprintf(
      'cat %s >> %s',
      $source_dir.'/info.txt_'.$_,
      $target_dir.'/info.txt'
    )
  ) for ('variation_converted', grep {$self->param($_)} qw(regulation));

  # create tar
  $self->tar(
    $self->dump_dir,
    $species,
    $self->param('assembly'),
    $suffix,
    '_tabixconverted'
  );
}

sub tar {
  my ($self, $dir, $species, $assembly, $suffix, $mod) = @_;

  $mod ||= '';

  my $tar_file = $self->get_tar_file_name($dir, $species, $assembly, $mod);
  
  unlink($tar_file);

  $self->run_system_command(
    # use some steps to prevent tar failing with warnings about files changing as they are read
    # --warning=no-file-changed suppresses that warning
    # "|| [[ $? -eq 1]]" allows us to ignore the return code of 1 which is returned if those warnings appear
    sprintf(
      'tar --warning=no-file-changed -cz -C %s -f %s %s || [[ $? -eq 1 ]]',
      $dir,
      $tar_file,
      $suffix
    )
  );
}

sub copy_synonyms {
  my $self = shift;
  my $target_dir = shift;

  $self->link_file(
    sprintf(
      '%s/synonyms/%s_%s_chr_synonyms.txt',
      $self->param('pipeline_dump_dir'),
      $self->param('species'),
      $self->param('assembly')
    ),
    $target_dir.'/chr_synonyms.txt'
  );
}

sub generate_md5s {
  my $self = shift;
  my $data_dir = $self->data_dir;
  my $data_dir_sed  = $data_dir;
  $data_dir_sed =~ s/\//\\\//g;
  my $command = 'find ' . $data_dir . ' -type f -exec md5sum "{}" + | sed \'s/' . $data_dir_sed . '\///\' > ' . $data_dir . '/MD5SUMS';
  $self->run_system_command($command);
}

1;
