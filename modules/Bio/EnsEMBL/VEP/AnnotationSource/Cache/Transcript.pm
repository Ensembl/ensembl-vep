=head1 LICENSE

Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript - local disk transcript annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseSerialized
  Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    gencode_basic
    all_refseq
    sift
    polyphen
    everything
    use_transcript_ref
    nearest
  )]);

  $self->check_sift_polyphen();

  # generate tree here otherwise forked processes will regenerate
  if($self->{nearest}) {
    throw("ERROR: --nearest requires Set::IntervalTree perl module to be installed\n") unless $Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript::CAN_USE_INTERVAL_TREE;
    $self->transcript_tree();
  }

  return $self;
}

sub check_sift_polyphen {
  my $self = shift;

  my $info = $self->info();

  foreach my $tool(qw(SIFT PolyPhen)) {
    my $lc_tool = lc($tool);
    my $v = $info->{$lc_tool};

    if($self->{$lc_tool}) {

      unless($v) {

        # dont die if user set "everything" param on a species with no SIFT/PolyPhen
        if($self->{everything} || $self->param('safe')) {
          $self->status_msg("INFO: disabling $tool");
          $self->param($lc_tool, 0);
          $self->{$lc_tool} = 0;
        }

        else {
          throw("ERROR: $tool not available\n");
        }
      }
    }
  }

  return 1;
}

sub get_dump_file_name {
  my $self = shift;
  my $chr  = shift;
  my $region = shift;

  throw("No chromosome given") unless $chr;
  throw("No region given") unless $region;

  # allow to pass region (start-end) or $start, $end
  $region .= '-'.shift if @_;

  return sprintf(
    "%s/%s/%s\.%s",
    $self->dir,
    $chr,
    $region,
    $self->file_suffix
  );
}

sub deserialized_obj_to_features {
  my $self = shift;
  my $obj = shift;

  my $tra = $self->get_adaptor('core', 'Translation');
  my $sa  = $self->get_adaptor('core', 'Slice');

  my @features;

  foreach my $chr(keys %$obj) {
    foreach my $t(@{$obj->{$chr}}) {

      # we may need to filter some out based on config options
      # do this ASAP to avoid processing more features than we need to
      next unless $self->filter_transcript($t);

      if(defined($t->{translation})) {
        $t->{translation}->{adaptor} = $tra;
        $t->{translation}->{transcript} = $t;
        weaken($t->{translation}->{transcript});
      }

      $t->{slice}->{adaptor} = $sa;

      $_->{slice} ||= $t->{slice} for @{$t->{_trans_exon_array}};

      push @features, $t;
    }
  }

  return \@features;
}

sub tree_file {
  my ($self, $tree) = @_;

  my $file = $self->_tree_coords_filename;

  # generate it by scanning the cache
  unless(-e $file) {
    my $as_dir = $self->dir;
    my $cache_region_size = $self->{cache_region_size};

    $self->status_msg("INFO: Scanning cache for transcript coordinates. This may take a while but will only run once per cache.\n");

    open TR, ">$file" or throw("ERROR: Could not write to transcript coords file: $!");

    opendir DIR, $as_dir;
    foreach my $c(grep {!/^\./ && -d $as_dir.'/'.$_} readdir DIR) {
      
      opendir CHR, $as_dir.'/'.$c;

      # scan each chr directory for transcript cache files e.g. 1-1000000.gz
      foreach my $file(grep {/\d+\-\d+\.gz/} readdir CHR) {
        my ($s) = split(/\D/, $file);

        # we can use parent classes methods to fetch transcripts into memory
        foreach my $t(
          grep {$_->biotype eq 'protein_coding'}
          @{$self->get_features_by_regions_uncached([[$c, ($s - 1) / $cache_region_size]])}
        ) {
          print TR "$c\t".join("\t", @{$self->_tree_file_data($t)})."\n";
        }

        # calling get_features_by_regions_uncached adds data to an in-memory cache
        # we need to clear it manually here
        $self->clean_cache;
      }
      closedir CHR;
    }
    closedir DIR;
    close TR;
  }

  return $file;
}

1;