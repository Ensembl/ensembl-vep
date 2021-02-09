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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::FinishDump;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::BaseVEP);

use File::Path qw(rmtree);

sub param_defaults {
  return {
    'dir_suffix' => '',
  };
}

sub run {
  my $self = shift;
  my $type = $self->param('type');

  if($type eq 'core') {
    $type = '';
  }
  else {
    $type = '_'.$type;
  }

  $self->rm_dirs($type);
  return;
}

sub rm_dirs {
  my $self = shift;
  my $type = shift;
  
  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');
  my $version  = $self->param('eg_version') || $self->required_param('ensembl_release');
  my $dir      = $self->data_dir;

  rmtree(
    sprintf(
      '%s/%s%s/%s_%s',
      $dir,
      $species,
      $type,
      $version,
      $assembly
    )
  );
  
  # remove species dir if empty
  my $species_dir = "$dir/$species$type";
  
  if(opendir DIR, $species_dir) {
    my @list = grep {!/^\.+$/} readdir DIR;
    closedir DIR;
    
    rmtree($species_dir) unless scalar @list;
  }

  # and the same in the dumps/ dir
  $dir = $self->dump_dir;
  $species_dir = "$dir/$species$type";
  
  if(opendir DIR, $species_dir) {
    my @list = grep {!/^\.+$/} readdir DIR;
    closedir DIR;
    
    rmtree($species_dir) unless scalar @list;
  }
  
  return;
}

1;
