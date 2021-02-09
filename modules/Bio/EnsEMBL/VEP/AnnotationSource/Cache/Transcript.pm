=head1 LICENSE

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript - local disk transcript annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
  config => $config,
  dir    => $dir
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Cache-based annotation source for transcript data.

Data are stored as serialized objects on disk. Structure:

$cache = {
  chr => [
    tr1,
    tr2,
    tr3
  ]
}

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript;

use Scalar::Util qw(weaken);
use File::Copy;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseSerialized
  Bio::EnsEMBL::VEP::AnnotationType::Transcript
);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config $config,
                 dir    => string $dir,
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript
  Exceptions : throws if "nearest" param set and Set::IntervalTree not installed
  Caller     : CacheDir
  Status     : Stable

=cut

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
    throw("ERROR: --nearest requires Set::IntervalTree perl module to be installed\n") unless $Bio::EnsEMBL::VEP::AnnotationType::Transcript::CAN_USE_INTERVAL_TREE;
    $self->transcript_tree();
  }

  return $self;
}


=head2 check_sift_polyphen

  Example    : $ok = $as->check_sift_polyphen();
  Description: Gets user-set SIFT/PolyPhen parameters and checks vs
               availability in the cache. If using "safe" mode (REST, web)
               or --everything, params are disabled if unavailable. Otherwise,
               this method will throw.
  Returntype : bool
  Exceptions : throws if configured tool not available
  Caller     : new()
  Status     : Stable

=cut

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


=head2 get_dump_file_name

  Arg 1      : string $chr
  Arg 2      : string $region or int $start
  Arg 3      : (optional) int $end
  Example    : $file = $as->get_dump_file_name(1, "1-1000000");
  Description: Gets file name from the cache given a region.
  Returntype : string
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

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


=head2 deserialized_obj_to_features

  Arg 1      : hashref $obj
  Example    : $features = $as->deserialized_obj_to_features($obj);
  Description: Takes the deserialized object read from the cache file,
               restores the objects by reattaching necessary adaptors.
  Returntype : arrayref of Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

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


=head2 tree_file

  Example    : $file = $as->tree_file();
  Description: Gets the filename storing data for populating the transcript
               tree for this annotation source. If it does not exist, then the
               file is created by scanning the cache for transcripts (this can
               take a while and a status message is printed indicating this).
  Returntype : string
  Exceptions : none
  Caller     : populate_tree()
  Status     : Stable

=cut

sub tree_file {
  my $self = shift;

  my $file = $self->_tree_coords_filename;

  # generate it by scanning the cache
  unless(-e $file) {
    my $as_dir = $self->dir;
    my $cache_region_size = $self->{cache_region_size};

    $self->status_msg("INFO: Scanning cache for transcript coordinates. This may take a while but will only run once per cache.\n");

    # use an intermediate tmp file so that if write is interrupted next run will start from scratch
    my $tmpfile = $file.".tmp".$$;
    open TR, ">$tmpfile" or throw("ERROR: Could not write to transcript coords file: $!");

    opendir DIR, $as_dir;
    foreach my $c(sort { $a cmp $b } grep {!/^\./ && -d $as_dir.'/'.$_} readdir DIR) {
      
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

    move($tmpfile, $file);
  }

  return $file;
}

1;
