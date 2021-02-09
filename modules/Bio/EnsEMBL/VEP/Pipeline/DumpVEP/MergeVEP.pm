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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::MergeVEP;

use strict;
use warnings;

use File::Path qw(make_path rmtree);
use Storable qw(nstore_fd fd_retrieve);

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

  my $data_dir = $self->data_dir;

  # do merging
  my $ens_root = $data_dir.'/'.$self->species_suffix;
  my $ref_root = $data_dir.'/'.$self->refseq_species_suffix;
  my $mrg_root = $data_dir.'/'.$self->merged_species_suffix;
  
  # clear any previous existing dir
  rmtree($mrg_root);

  opendir ENSROOT, $ens_root;
  opendir REFROOT, $ref_root;

  # list chromosomes
  my %ens_chrs = map {$_ => 1} grep {-d $ens_root.'/'.$_ && !/^\./} readdir ENSROOT;
  my %ref_chrs = map {$_ => 1} grep {-d $ref_root.'/'.$_ && !/^\./} readdir REFROOT;

  closedir ENSROOT;
  closedir REFROOT;

  # mkdirs
  my %mrg_chrs = map {$_ => 1} (keys %ens_chrs, keys %ref_chrs);

  foreach my $chr(keys %mrg_chrs) {
    make_path($mrg_root.'/'.$chr);
  
    # exists in both
    if(-d $ens_root.'/'.$chr && -d $ref_root.'/'.$chr) {
      opendir ENSCHR, $ens_root.'/'.$chr;
      opendir REFCHR, $ref_root.'/'.$chr;
    
      my %ens_files = map {$_ => 1} grep {/\d+\.gz$/} readdir ENSCHR;
      my %ref_files = map {$_ => 1} grep {/\d+\.gz$/} readdir REFCHR;
    
      closedir ENSCHR;
      closedir REFCHR;
    
      my %mrg_files = map {$_ => 1} (keys %ens_files, keys %ref_files);
    
      foreach my $file(keys %mrg_files) {
      
        # exists in both, need to concatenate
        if(-e $ens_root.'/'.$chr.'/'.$file && -e $ref_root.'/'.$chr.'/'.$file) {

          # read in Ensembl cache
          open my $ens_fh, "gzip -dc ".$ens_root.'/'.$chr.'/'.$file." |";
          my $ens_cache;
          $ens_cache = fd_retrieve($ens_fh);
          close $ens_fh;
        
          # add a flag to each transcript indicating which cache it came from
          $_->{_source_cache} = 'Ensembl' for @{$ens_cache->{$chr}};
        
          # do same for RefSeq
          open my $ref_fh, "gzip -dc ".$ref_root.'/'.$chr.'/'.$file." |";
          my $ref_cache;
          $ref_cache = fd_retrieve($ref_fh);
          close $ref_fh;
          $_->{_source_cache} = 'RefSeq' for @{$ref_cache->{$chr}};
        
          # merge and sort transcript lists
          my $mrg_cache;
          @{$mrg_cache->{$chr}} = sort {$a->{start} <=> $b->{start}} (@{$ens_cache->{$chr}}, @{$ref_cache->{$chr}});
        
          # dump to new file
          open my $mrg_fh, "| gzip -9 -c > ".$mrg_root.'/'.$chr.'/'.$file or die "ERROR: Could not write to dump file";
          nstore_fd($mrg_cache, $mrg_fh);
          close $mrg_fh;
        }
      
        # otherwise simply copy/link
        else {
          my $root = -e $ens_root.'/'.$chr.'/'.$file ? $ens_root : $ref_root;
          $self->link_file($root.'/'.$chr.'/'.$file, $mrg_root.'/'.$chr.'/'.$file);
        }
      }
    }
  
    # only exists in one, simply copy all files
    else {
      my $root = -d $ens_root.'/'.$chr ? $ens_root : $ref_root;
    
      opendir CHR, $root.'/'.$chr;
      $self->link_file($root.'/'.$chr.'/'.$_, $mrg_root.'/'.$chr.'/'.$_) for grep {!/^\./} readdir CHR;
      closedir CHR;
    }
  }

  unlink($mrg_root.'/info.txt_core');
  $self->run_system_command(
    sprintf(
      'cat %s %s | sort -u >> %s',
      $ens_root.'/info.txt_core',
      $ref_root.'/info.txt_core',
      $mrg_root.'/info.txt_core'
    )
  );
}


1;
