=head1 LICENSE

Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSourceAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSourceAdaptor - gets all AnnotationSources from initial config

=head1 SYNOPSIS

my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({
  config => $config
});

$sources = $asa->get_all();

=head1 DESCRIPTION

Factory for generating AnnotationSources from configuration in supplied Bio::EnsEMBL::VEP::Config.

Can create database- and file-based sources directly; uses Bio::EnsEMBL::VEP::CacheDir to generate
cache-based sources.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSourceAdaptor;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::CacheDir;
use Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat;
use Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation;
use Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation;
use Bio::EnsEMBL::VEP::AnnotationSource::File;

use LWP::Simple;

=head2 get_all

  Example    : $sources = $asa->get_all()
  Description: Gets all AnnotationSources
  Returntype : arrayref of Bio::EnsEMBL::VEP::AnnotationSource
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::BaseRunner
  Status     : Stable

=cut

sub get_all {
  my $self = shift;

  return [
    sort {($b->{can_filter_vfs} || 0) <=> ($a->{can_filter_vfs} || 0)}
    (
      @{$self->get_all_from_cache},
      @{$self->get_all_from_database},
      @{$self->get_all_custom},
    )
  ];
}


=head2 get_all_from_cache

  Example    : $sources = $asa->get_all_from_cache()
  Description: Gets all cache AnnotationSources
  Returntype : arrayref of Bio::EnsEMBL::VEP::AnnotationSource
  Exceptions : none
  Caller     : get_all()
  Status     : Stable

=cut

sub get_all_from_cache {
  my $self = shift;

  return [] unless $self->param('cache');

  my $cache_dir_obj = Bio::EnsEMBL::VEP::CacheDir->new({
    config   => $self->config,
    root_dir => $self->param('dir_cache') || $self->param('dir'),
    dir      => $self->param('full_cache_dir')
  });

  return $cache_dir_obj->get_all_AnnotationSources();
}


=head2 get_all_from_database

  Example    : $sources = $asa->get_all_from_database()
  Description: Gets all database AnnotationSources
  Returntype : arrayref of Bio::EnsEMBL::VEP::AnnotationSource
  Exceptions : none
  Caller     : get_all()
  Status     : Stable

=cut

sub get_all_from_database {
  my $self = shift;

  return [] if $self->param('offline');

  my @as;

  # we don't want to get e.g. transcript DB sources if we have cache
  unless($self->param('cache') || ($self->param('custom') && !$self->param('database'))) {
    my $module = $self->module_prefix.'::AnnotationSource::Database::Transcript';
    eval "require $module";

    if($self->param('database')) {
      push @as, $module->new({
        config => $self->config,
        filter => $self->param('transcript_filter'),
        bam    => $self->param('bam'),
      });

      # special case merged
      if($self->param('merged') && $self->get_adaptor('otherfeatures', 'slice')) {
        my $core_type_bak = $self->param('core_type');
        $self->param('core_type', 'otherfeatures');

        push @as, $module->new({
          config => $self->config,
          filter => $self->param('transcript_filter'),
          bam    => $self->param('bam'),
        });

        $self->param('core_type', $core_type_bak);
      }
    }

    push @as, Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat->new({
      config => $self->config,
    }) if $self->param('regulatory') && $self->get_adaptor('funcgen', 'RegulatoryFeature');

    push @as, Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
      config => $self->config,
    }) if $self->param('check_existing') && $self->get_adaptor('variation', 'Variation');
  }

  # overlapping SVs
  # this has no cache equivalent
  push @as, Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation->new({
    config => $self->config,
  }) if $self->param('check_svs') && $self->get_adaptor('variation', 'Variation');

  return \@as;
}


=head2 get_all_custom

  Example    : $sources = $asa->get_all_custom()
  Description: Gets all custom file AnnotationSources
  Returntype : arrayref of Bio::EnsEMBL::VEP::AnnotationSource
  Exceptions : none
  Caller     : get_all()
  Status     : Stable

=cut

sub get_all_custom {
  my $self = shift;

  my @as;

  my @VALID_OPTIONS = (
    'file',
    'format',
    'short_name',
    'fields',
    'type',
    'overlap_cutoff',
    'reciprocal',
    'distance',
    'coords',
    'same_type',
    'num_records',
    'summary_stats'
  );

  foreach my $custom_string(@{$self->param('custom') || []}) {
    
    my %hash = ();
    my (@params, @fields);

    if (grep(/=/, $custom_string)){
      @params = split /\,/, $custom_string;
    } else {
      ($hash{"file"}, $hash{"short_name"}, $hash{"format"}, $hash{"type"}, $hash{"coords"}, @fields) = split /\,/, $custom_string;
    }

    if (@params){
      foreach my $param(@params) {
        my ($key, $val) = split('=', $param);
        die("ERROR: Failed to parse parameter $param; Please add <VALUE_OF_PARAMETER>=$param\n") unless defined($key) && defined($val);
        $hash{$key} = $val;
      };

      # warn about invalid options
      my @invalid_opts;
      for my $opt (keys %hash) {
        push @invalid_opts, $opt unless grep { $opt eq $_ } @VALID_OPTIONS;
      }
      throw("ERROR: The following options are not supported for custom annotations: "
            . join(", ", @invalid_opts) . "\n" .
            "LINE: --custom $custom_string\n") if @invalid_opts;
    };

    throw("ERROR: No file was added for custom annotation source.\nLINE: --custom $custom_string\n") unless defined($hash{"file"});
    throw("ERROR: No format specified for custom annotation source.\nLINE: --custom $custom_string\n") unless defined($hash{"format"});
    throw("ERROR: Access to remote data files disabled\n") if $self->param('no_remote') && $hash{"file"} =~ /^(ht|f)tp:\/\/.+/;

    my $opts = {
      config => $self->config,
      file => $hash{"file"},
      short_name => $hash{"short_name"},
      format => $hash{"format"},
      type => $hash{"type"} || "overlap",
      report_coords => $hash{"coords"} || 0,
      overlap_cutoff => $hash{"overlap_cutoff"} || 0,
      distance => $hash{"distance"},
      same_type => $hash{"same_type"} || 0,
      reciprocal => $hash{"reciprocal"} || 0,
      num_records => $hash{"num_records"}
    };

    $opts->{overlap_def} = $opts->{reciprocal} ?
      "Percentage of minimum reciprocal overlap between input variant and reference variant" :
      "Percentage of input variant covered by reference variant";

    my $format = $hash{"format"};
    if(defined($format) && $format =~ /^G[TF]F$/i) {
      $opts->{filter} = $self->param('transcript_filter');
      $opts->{bam} = $self->param('bam');
    }

    $opts->{fields} = [split /%/, $hash{"fields"}] if $hash{"fields"};
    $opts->{fields} = \@fields if @fields;

    if (!defined $opts->{num_records}) {
      $opts->{num_records} = 50; # by default, show all values for non-SVs
    } elsif ($opts->{num_records} eq 'all') {
      $opts->{num_records} = 'inf';
    }

    # Default summary statistics: only show for BED/bigwig custom files
    $opts->{summary_stats} = $hash{"summary_stats"} || 'none';
    delete $opts->{summary_stats} if $opts->{summary_stats} eq 'none';

    if ( $opts->{summary_stats} ) {
      $opts->{summary_stats} = [split /%/, $opts->{summary_stats}];

      # Check invalid summary statistics
      my @invalid;
      my @stats = ('min', 'max', 'mean', 'count', 'sum');
      for my $k (@{ $opts->{summary_stats} }) {
        push @invalid, $k unless grep { $_ eq $k } @stats;
      }

      if (@invalid) {
        my $invalid_opts = join(", ", @invalid);
        my $valid_opts   = join(", ", @stats);
        throw("ERROR: The following summary statistics for custom annotations ".
              "are not supported: $invalid_opts. ".
              "Available options are: $valid_opts.\n".
              "LINE: --custom $custom_string\n");
      }
    }

    if (grep { /\#\#\#CHR\#\#\#/ } $hash{"file"}){

      my @valid_chromosomes = keys %{$self->chr_lengths} > 0 ? sort keys %{$self->chr_lengths}: ((1..22), qw(X Y MT));
      
      foreach my $chr (@valid_chromosomes){
        print $chr."\n";
        my $new_file = $hash{"file"};
        my $new_opts = { %$opts };
        $new_file =~ s/\#\#\#CHR\#\#\#/$chr/;
        next unless ( -e $new_file || head($new_file) );
        $new_opts->{file} = $new_file;
        push @as, Bio::EnsEMBL::VEP::AnnotationSource::File->new($new_opts);
      }

      # Non-match ###CHR### pattern scenario
      die "Error: No files with pattern " . $hash{"file"} . " were found\n" unless @as;
    
    } else {
      push @as, Bio::EnsEMBL::VEP::AnnotationSource::File->new($opts);
    }
  }

  return \@as;
}

1;
