=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::CacheDir;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix;

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

# this method is called from VEP::Runner's init() method
sub get_all_AnnotationSources {
  my $self = shift;

  if(!exists($self->{_annotation_sources})) {

    my $dir = $self->dir;
    my $info = $self->info;

    my @as;

    # initialise with transcript source
    push @as, Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
      config => $self->config,
      dir => $dir,
      serializer_type => $info->{serialiser_type} || undef,
      source_type => $self->source_type,
      cache_region_size => $info->{cache_region_size} || $self->param('cache_region_size'),
      info => $self->version_data,
    });

    # add RegFeats if available
    push @as, Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat->new({
      config => $self->config,
      dir => $dir,
      serializer_type => $info->{serialiser_type} || undef,
      cache_region_size => $info->{cache_region_size} || $self->param('cache_region_size'),
      info => $self->version_data,
    }) if $self->param('regulatory') and $info->{regulatory};

    # add Variation if available
    if($self->param('check_existing') && $info->{variation_cols}) {

      my $class = 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation';

      # tabix type?
      if(($info->{var_type} || '') eq 'tabix') {
        $class .= 'Tabix';
      }

      push @as, $class->new({
        config => $self->config,
        dir => $dir,
        cache_region_size => $info->{cache_region_size} || $self->param('cache_region_size'),
        cols => $info->{variation_cols},
        info => $self->version_data,
      }); 
    }

    $self->{_annotation_sources} = \@as;
  }

  return $self->{_annotation_sources};
}

sub init {
  my $self = shift;

  # check if there's a FASTA file in there
  unless($self->param('fasta') || $self->param('no_fasta')) {
    my $dir = $self->dir;

    opendir CACHE, $dir;
    my ($fa) = grep {/\.fa(\.gz)?$/} readdir CACHE;
    closedir CACHE;

    if(defined $fa) {
      $self->param('fasta', $dir.'/'.$fa);
      $self->status_msg("Auto-detected FASTA file in cache directory") if $self->param('verbose');
    }
  }

  return 1;
}

sub root_dir {
  my $self = shift;
  $self->{root_dir} = shift if @_;
  return $self->{root_dir};
}

sub dir {
  my $self = shift;
  
  $self->{dir} = shift if @_;

  if(!exists($self->{dir})) {

    # initialise dir at root_dir
    my $dir = $self->root_dir;

    # complete dir with species name and db_version
    my $species_dir_name = $self->species();
    $species_dir_name .= '_'.$_ for grep { $self->param($_) } qw(refseq merged);

    # add species dir name
    $dir .= '/'.$species_dir_name;

    # check whats in here to match to assembly if given
    die("ERROR: Cache directory $dir not found\n") if !-e $dir;

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
            (
              $self->param('host') eq 'ensembldb.ensembl.org' ?
              "\nIf using human GRCh37 add \"--port 3337\"".
              " to use the GRCh37 database, or --offline to avoid database connection entirely\n" :
              ''
            )
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

    ## NOT CURRENTLY BEING USED, COMMENTED OUT AS NO UNIT TEST
    # check if any disabled options are in use
    # these are set in the cache info file
    # if(my $cache_disabled = $info->{cache_disabled}) {
    #   my @arr = ref($cache_disabled eq 'ARRAY') ? @{$cache_disabled} : ($cache_disabled);

    #   if(my ($disabled) = grep {$self->param($_)} @arr) {
    #     throw("ERROR: Unable to use --".$disabled." with this cache\n");
    #   }
    # }

    $self->{info} = $info;
  }

  return $self->{info};
}

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

sub version_data {
  return $_[0]->info->{version_data} || {};
}

sub source_type {
  my $self = shift;

  if(!exists($self->{source_type})) {
    ($self->{source_type}) = grep {$self->param($_)} qw(merged refseq);
    $self->{source_type} ||= 'ensembl';
  }

  return $self->{source_type};
}

1;