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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSourceAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSourceAdaptor - gets all AnnotationSources from initial config

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

# this method is called from VEP::Runner's init() method
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

sub get_all_from_database {
  my $self = shift;

  return [] if $self->param('offline');

  my @as;

  # we don't want to get e.g. transcript DB sources if we have cache
  unless($self->param('cache') || ($self->param('custom') && !$self->param('database'))) {
    my $module = $self->module_prefix.'::AnnotationSource::Database::Transcript';
    eval "require $module";
    push @as, $module->new({
      config => $self->config,
      filter => $self->param('transcript_filter'),
    }) if $self->param('database');

    push @as, Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat->new({
      config => $self->config,
    }) if $self->param('regulatory');

    push @as, Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
      config => $self->config,
    }) if $self->param('check_existing');
  }

  # overlapping SVs
  # this has no cache equivalent
  push @as, Bio::EnsEMBL::VEP::AnnotationSource::Database::StructuralVariation->new({
    config => $self->config,
  }) if $self->param('check_svs');

  return \@as;
}

sub get_all_custom {
  my $self = shift;

  my @as;

  foreach my $custom_string(@{$self->param('custom') || []}) {
    my ($file, $short_name, $format, $type, $report_coords, @fields) = split /\,/, $custom_string;

    throw("ERROR: No format specified for custom annotation source $file\n") unless $format;

    my $opts = {
      config => $self->config,
      file => $file,
      short_name => $short_name,
      format => $format,
      type => $type,
      report_coords => $report_coords,
    };

    $opts->{filter} = $self->param('transcript_filter') if $format =~ /^G[TF]F$/i;

    $opts->{fields} = \@fields if @fields;
    
    push @as, Bio::EnsEMBL::VEP::AnnotationSource::File->new($opts);
  }

  return \@as;
}

1;