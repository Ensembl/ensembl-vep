=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript - local disk transcript annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::AnnotationSource::Cache::Transcript;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::Haplo::AnnotationSource::BaseTranscript Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript);

sub populate_tree {
  my ($self, $tree) = @_;

  my $as_dir = $self->dir;

  if(-e $as_dir.'/transcript_coords.txt') {
    open TR, $as_dir.'/transcript_coords.txt';
    while(<TR>) {
      chomp;
      $tree->insert(split);
    }
    close TR;
  }
  else {
    my $cache_region_size = $self->{cache_region_size};

    open TR, ">".$as_dir.'/transcript_coords.txt' or throw("ERROR: Could not write to transcript coords file: $!");

    opendir DIR, $as_dir;
    foreach my $c(grep {!/^\./ && -d $as_dir.'/'.$_} readdir DIR) {
      
      opendir CHR, $as_dir.'/'.$c;
      foreach my $file(grep {/\d+\-\d+\.gz/} readdir CHR) {
        my ($s) = split(/\D/, $file);

        foreach my $t(
          grep {$_->biotype eq 'protein_coding'}
          @{$self->get_features_by_regions_uncached([[$c, ($s - 1) / $cache_region_size]])}
        ) {
          my ($s, $e) = ($t->seq_region_start, $t->seq_region_end);
          $tree->insert($c, $s, $e);
          print TR "$c\t$s\t$e\n";
        }
      }
      closedir CHR;
    }
    closedir DIR;
    close TR;
  }
}

1;