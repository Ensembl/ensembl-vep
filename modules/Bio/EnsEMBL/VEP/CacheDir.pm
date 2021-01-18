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

# EnsEMBL module for Bio::EnsEMBL::VEP::CacheDir
#
#

=head1 NAME

Bio::EnsEMBL::VEP::CacheDir - a class representing a "traditional" VEP cache directory

=head1 SYNOPSIS

my $cache_dir = Bio::EnsEMBL::VEP::CacheDir->new({
  config   => $config,
  root_dir => '/nfs/home/joebloggs/.vep/'
});

my $info = $cache_dir->info();

my $sources = $cache_dir->get_all_AnnotationSources()

=head1 DESCRIPTION

The CacheDir class represents a VEP cache directory containing
transcript data and potentially also regulatory and variation
data.

The class contains methods for determining which data is available,
reading metadata, and generating AnnotationSource objects.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::CacheDir;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix;
use Bio::EnsEMBL::Variation::Utils::FastaSequence;


=head2 new

  Arg 1      : hashref $args
               {
                 config   => Bio::EnsEMBL::VEP::Config,
                 dir      => string (e.g. /nfs/home/joebloggs/.vep/homo_sapiens/88_GRCh38),
                 root_dir => string (e.g. /nfs/home/joebloggs/.vep/),
               }
  Example    : $obj = Bio::EnsEMBL::VEP::CacheDir->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::CacheDir object. If a root_dir
               arg is supplied, the full dir is filled out from various
               config parameters.
  Returntype : Bio::EnsEMBL::VEP::CacheDir
  Exceptions : throws if neither dir or root_dir given
  Caller     : Bio::EnsEMBL::VEP::AnnotationSourceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  # need to either specify full path (dir) e.g. ~/.vep/homo_sapiens/78_GRCh38/
  # or root path (root_dir) e.g. ~/.vep/ where species and version are drawn from config and assembly from config, DB or scanning dir
  throw("ERROR: No root_dir or dir specified") unless $hashref && ($hashref->{root_dir} or $hashref->{dir});

  $self->{$_} = $hashref->{$_} for keys %$hashref;

  $self->init();

  return $self;
}


=head2 get_all_AnnotationSources

  Example    : $sources = $cache_dir->get_all_AnnotationSources();
  Description: Get all AnnotationSources from this CacheDir given
               user parameters from config.
  Returntype : listref of Bio::EnsEMBL::VEP::AnnotationSource
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::AnnotationSourceAdaptor
  Status     : Stable

=cut

sub get_all_AnnotationSources {
  my $self = shift;

  if(!exists($self->{_annotation_sources})) {

    my $dir = $self->dir;
    my $info = $self->info;

    my @as;

    # initialise with transcript source
    my $module = $self->module_prefix.'::AnnotationSource::Cache::Transcript';
    eval "require $module";
    push @as, $module->new({
      config => $self->config,
      dir => $dir,
      serializer_type => $info->{serialiser_type} || undef,
      source_type => $self->source_type,
      cache_region_size => $info->{cache_region_size} || $self->param('cache_region_size'),
      info => $self->version_data,
      valid_chromosomes => $info->{valid_chromosomes},
      filter => $self->param('transcript_filter'),
      bam => $self->param('bam'),
    });

    # add RegFeats if available
    push @as, Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat->new({
      config => $self->config,
      dir => $dir,
      serializer_type => $info->{serialiser_type} || undef,
      cache_region_size => $info->{cache_region_size} || $self->param('cache_region_size'),
      info => $self->version_data,
      available_cell_types => [split(',', ($info->{cell_types} || ''))],
      valid_chromosomes => $info->{valid_chromosomes},
    }) if $self->param('regulatory') and $info->{regulatory};

    # add Variation if available
    if($self->param('check_existing') && $info->{variation_cols}) {

      my $class = 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation';

      # tabix type?
      if(($info->{var_type} || '') eq 'tabix') {
        $class .= 'Tabix';
      }

      # check for presence of ExAC/gnomAD vs what user has requested
      # don't do this check with --everything
      unless($self->param('everything')) {
        my $have_exac   = (grep {$_ eq 'ExAC'} @{$info->{variation_cols}});
        my $have_gnomad = (grep {$_ eq 'gnomAD'} @{$info->{variation_cols}});

        if($self->param('af_exac') && !$have_exac) {
          throw(
            "ERROR: ExAC data is not available in this cache".
            ($have_gnomad ? "; gnomAD exome data is available with --af_gnomad\n" : "\n")
          );
        }
        if($self->param('af_gnomad') && !$have_gnomad) {
          throw(
            "ERROR: gnomAD data is not available in this cache".
            ($have_exac ? "; ExAC data is available with --af_exac\n" : "\n")
          );
        }
      }

      push @as, $class->new({
        config => $self->config,
        dir => $dir,
        cache_region_size => $info->{cache_region_size} || $self->param('cache_region_size'),
        cols => $info->{variation_cols},
        info => $self->version_data,
        valid_chromosomes => $info->{valid_chromosomes},
      }); 
    }

    $self->{_annotation_sources} = \@as;
  }

  return $self->{_annotation_sources};
}


=head2 init

  Example    : $cache_dir->init();
  Description: Initialise the CacheDir object: complete the full dir path,
               auto-detect FASTA and chromosome synonyms files.
  Returntype : bool
  Exceptions : none
  Caller     : new()
  Status     : Stable

=cut

sub init {
  my $self = shift;

  my $dir = $self->dir;

  # check if there's a FASTA file in there
  unless($self->param('fasta') || $self->param('no_fasta')) {
    opendir CACHE, $dir;

    # look for a .fa or a .fa.gz depending on whether we have Bio::DB::HTS installed
    my $fa;
    if($Bio::EnsEMBL::Variation::Utils::FastaSequence::CAN_USE_FAIDX) {
      ($fa) = grep {/\.fa(\.gz)?$/} readdir CACHE;
    }
    else {
      ($fa) = grep {/\.fa$/} readdir CACHE; 
    }
    closedir CACHE;

    if(defined $fa) {
      $self->param('fasta', $dir.'/'.$fa);
      $self->status_msg("Auto-detected FASTA file in cache directory") if $self->param('verbose');
    }
  }

  $self->param('synonyms', $dir.'/chr_synonyms.txt') if !$self->param('synonyms') && -e $dir.'/chr_synonyms.txt';

  return 1;
}


=head2 root_dir
  
  Arg 1      : (optional) string $root_dir
  Example    : $root_dir = $cache_dir->root_dir();
  Description: Getter/setter for root_dir field; this is the root directory
               in which a user's VEP caches are installed.
  Returntype : string
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub root_dir {
  my $self = shift;
  $self->{root_dir} = shift if @_;
  return $self->{root_dir};
}


=head2 dir
  
  Arg 1      : (optional) string $dir
  Example    : $dir = $cache_dir->dir();
  Description: Getter/setter for dir field; this is the full path for a specific
               species/version/assembly cache.
  Returntype : string
  Exceptions : throws if:
                - [root_dir]/[species] directory does not exist
                - [species]/[version]_* does not exist
                - [species]/[version]_[assembly] exists but [assembly] does not match assembly param
                - multiple assemblies found for [version] but none specified or inferred from database
                - [species]/[version]_[assembly] does not exist
  Caller     : internal
  Status     : Stable

=cut

sub dir {
  my $self = shift;
  
  $self->{dir} = shift if @_;

  if(!defined($self->{dir})) {

    # initialise dir at root_dir
    my $dir = $self->root_dir;

    # complete dir with species name and db_version
    my $species_dir_name = $self->species();
    $species_dir_name .= '_'.$_ for grep { $self->param($_) } qw(refseq merged);

    # add species dir name
    $dir .= '/'.$species_dir_name;

    # check whats in here to match to assembly if given
    throw("ERROR: Cache directory $dir not found\n") if !-e $dir;

    # get version
    # user may specify this with --cache_version or --db_version
    # if they don't, use registry software version (maybe change this to VEP version?)
    my $cache_version = $self->param('cache_version') || $self->param('db_version') || $self->registry->software_version;

    opendir DIR, $dir;
    my @dir_contents = grep {!/^\./} readdir DIR;
    closedir DIR;

    my @matched_contents = grep {/^$cache_version/} @dir_contents;

    # find out what the user specified as the assembly version
    my $config_assembly = $self->param('assembly');

    # no matched entries, cache not installed
    if(scalar @matched_contents == 0) {
      throw("ERROR: No cache found for $species_dir_name, version $cache_version\n");
    }

    # only 1 entry, can assume this is OK
    elsif(scalar @matched_contents == 1) {
      $dir .= '/'.$matched_contents[0];

      # is there an assembly here?
      if($matched_contents[0] =~ /\d+\_(.+)/) {
        my $matched_assembly = $1;

        if(defined($config_assembly) && $config_assembly ne $matched_assembly) {
          throw(
            "ERROR: Cache assembly version ($matched_assembly) ".
            "and database or selected assembly version (".$config_assembly.
            ") do not match\n".
            "\n You can use --offline to avoid a database connection entirely\n" 
          );
        }

        $self->param('assembly', $matched_assembly) unless $config_assembly;
      }
    }

    # did user specify assembly version?
    elsif(!defined($config_assembly)) {
      my $possibles = join(", ", map {s/^$cache_version\_//; $_} @matched_contents);
      throw("ERROR: Multiple assemblies found for cache version $cache_version ($possibles) - specify one using --assembly [assembly]\n");
    }

    # add cache version and assembly
    else {
      $dir .= '/'.$cache_version.'_'.$config_assembly;
      throw("ERROR: No cache found for ".$self->species.", version $cache_version, assembly ".$config_assembly."\n") unless -e $dir;
    }

    # nuke consecutive "/"
    $dir =~ s/\/+/\//g;

    $self->{dir} = $dir;
  }

  return $self->{dir};
}


=head2 info
  
  Arg 1      : (optional) hashref $info
  Example    : $info = $cache_dir->info();
  Description: Getter/setter for info hashref. If not set, reads metadata from
               cache info file
  Returntype : hashref
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub info {
  my $self = shift;
  $self->{info} = shift if @_;

  if(!exists($self->{info})) {

    # read cache info
    my $info = $self->read_cache_info_file($self->dir);

    my $config_assembly = $self->param('assembly');

    # check assembly matches
    if(defined($config_assembly) && defined($info->{assembly}) && $config_assembly ne $info->{assembly}) {
      throw("ERROR: Mismatch in assembly versions from config (".$config_assembly.") and cache info.txt file (".$info->{assembly}.")\n");
    }

    # enable use_transcript_ref
    if($info->{bam} && !$self->param('use_given_ref')) {
      $self->status_msg("INFO: BAM-edited cache detected, enabling --use_transcript_ref; use --use_given_ref to override this\n");
      $self->param('use_transcript_ref', 1);
      $self->param('bam_edited', 1);
    }

    ## NOT CURRENTLY BEING USED, COMMENTED OUT AS NO UNIT TEST
    # check if any disabled options are in use
    # these are set in the cache info file
    # if(my $cache_disabled = $info->{cache_disabled}) {
    #   my @arr = ref($cache_disabled eq 'ARRAY') ? @{$cache_disabled} : ($cache_disabled);

    #   if(my ($disabled) = grep {$self->param($_)} @arr) {
    #     throw("ERROR: Unable to use --".$disabled." with this cache\n");
    #   }
    # }

    # get valid chromosomes
    if(opendir DIR, $self->dir) {
      $info->{valid_chromosomes} = [sort grep {!/^\./ && -d $self->dir.'/'.$_} readdir DIR];
      closedir DIR;
    }

    $self->{info} = $info;
  }

  return $self->{info};
}


=head2 read_cache_info_file
  
  Arg 1      : string $dir
  Example    : $info = $cache_dir->read_cache_info_file($dir);
  Description: Reads cache metadata from [dir]/info.txt
  Returntype : hashref
  Exceptions : none
  Caller     : info()
  Status     : Stable

=cut

sub read_cache_info_file {
  my $self = shift;
  my $dir = shift;

  my $file = $dir.'/info.txt';

  my $cache_info = {};

  if(open IN, $file) {
    while(<IN>) {
      next if /^#/;
      chomp;
      my ($param, $value) = split "\t";

      if($param =~ s/^source_//) {
        $cache_info->{version_data}->{$param} = $value;
      }
      elsif($param =~ /variation_col/) {
        $cache_info->{$param} = [split ',', $value];
      }
      else {
        $cache_info->{$param} = $value unless defined $value && $value eq '-';
      }
    }

    close IN;
  }

  return $cache_info;
}


=head2 version_data
  
  Example    : $version_data = $cache_dir->version_data();
  Description: Gets version data for cache
  Returntype : hashref
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub version_data {
  return $_[0]->info->{version_data} || {};
}


=head2 source_type
  
  Example    : $source_type = $cache_dir->source_type();
  Description: Gets configured source type (ensembl, refseq or merged)
  Returntype : string
  Exceptions : none
  Caller     : dir()
  Status     : Stable

=cut

sub source_type {
  my $self = shift;

  if(!exists($self->{source_type})) {
    ($self->{source_type}) = grep {$self->param($_)} qw(merged refseq);
    $self->{source_type} ||= 'ensembl';
  }

  return $self->{source_type};
}

1;
