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

# EnsEMBL module for Bio::EnsEMBL::VEP::BaseRunner
#
#

=head1 NAME

Bio::EnsEMBL::VEP::BaseRunner - base runner class for VEP

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::BaseRunner;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use FileHandle;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(get_time merge_hashes get_version_data);
use Bio::EnsEMBL::VEP::Constants;
use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::InputBuffer;
use Bio::EnsEMBL::VEP::AnnotationSourceAdaptor;

# has our own new method, does not use BaseVEP's
# since this is the class users will be instantiating
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # initialise self
  my $self = bless {}, $class;

  # get a config object
  $self->{_config} = Bio::EnsEMBL::VEP::Config->new(@_);

  return $self;
}

sub setup_db_connection {
  my $self = shift;

  return if $self->param('offline');
  return unless $self->param('database') || $self->param('cache');

  # doing this inits the registry and DB connection
  my $reg = $self->registry();

  # check assembly
  if(my $db_assembly = $self->get_database_assembly) {

    my $config_assembly = $self->param('assembly');

    throw(
      "ERROR: Assembly version specified by --assembly (".$config_assembly.
      ") and assembly version in coord_system table (".$db_assembly.") do not match\n".
      (
        $self->param('host') eq 'ensembldb.ensembl.org' ?
        "\nIf using human GRCh37 add \"--port 3337\"".
        " to use the GRCh37 database, or --offline to avoid database connection entirely\n" :
        ''
      )
    ) if $config_assembly && $config_assembly ne $db_assembly;

    # update to database version
    $self->param('assembly', $db_assembly);

    if(!$self->param('assembly')) {
      throw("ERROR: No assembly version specified, use --assembly [version] or check the coord_system table in your core database\n");
    }
  }

  # update species, e.g. if user has input "human" we get "homo_sapiens"
  my $latin_species = $reg->get_alias($self->param('species'));
  # $latin_species =~ s/\d+$//;
  $self->species($latin_species);

  return 1;
}

sub get_all_AnnotationSources {
  my $self = shift;

  if(!exists($self->{_annotation_sources})) {
    my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $self->config});
    $self->{_annotation_sources} = $asa->get_all;
  }

  return $self->{_annotation_sources};
}

sub get_output_file_handle {
  my $self = shift;

  unless(exists($self->{output_file_handle})) {

    my $output_file_handle = FileHandle->new();

    my $output_file_name = $self->param('output_file');
      
    # check if file exists
    if(-e $output_file_name && !$self->param('force_overwrite')) {
      throw("ERROR: Output file $output_file_name already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }
    
    if(uc($output_file_name) eq 'STDOUT') {
      $output_file_handle = *STDOUT;
    }
    else {
      $output_file_handle->open(">$output_file_name") or throw("ERROR: Could not write to output file $output_file_name\n");
    }

    $self->{output_file_handle} = $output_file_handle;
  }

  return $self->{output_file_handle};
}

sub get_stats_file_handle {
  my $self = shift;
  my $type = shift;

  my $stats_file_root = $self->param('stats_file') || $self->param('output_file').'_summary';

  my $file_name;

  if($stats_file_root =~ m/^(.+?)\.(.+)$/) {
    my ($root, $ext) = ($1, $2);
    if(lc($ext) eq $type) {
      $file_name = $stats_file_root;
    }
    elsif(lc($ext) =~ /^(txt|html)$/) {
      $file_name = $root.'.'.$type;
    }
    else {
      $file_name = $stats_file_root.'.'.$type;
    }
  }
  else {
    $file_name = $stats_file_root.'.'.$type;
  }
      
  # check if file exists
  if(-e $file_name && !$self->param('force_overwrite')) {
    throw("ERROR: Stats file $file_name already exists. Specify a different output file with --stats_file or overwrite existing file with --force_overwrite\n");
  }
  
  my $fh = FileHandle->new();
  $fh->open(">$file_name") or throw("ERROR: Could not write to stats file $file_name\n");

  return $fh;
}

sub dump_stats {
  my $self = shift;

  unless($self->param('no_stats')) {

    if($self->param('stats_text') || !$Bio::EnsEMBL::VEP::Stats::CAN_USE_CGI) {
      my $fh = $self->get_stats_file_handle('txt');
      $self->stats->dump_text($fh);
      close $fh;
    }

    if($self->param('stats_html') && $Bio::EnsEMBL::VEP::Stats::CAN_USE_CGI) {
      my $fh = $self->get_stats_file_handle('html');
      $self->stats->dump_html($fh);
      close $fh;
    }
  }
}

sub get_output_header_info {
  my $self = shift;

  if(!exists($self->{output_header_info})) {

    my $vep_version_data = get_version_data()->{'ensembl-vep'};

    my $info = {
      time          => get_time,
      vep_version   => $vep_version_data->{release}.(defined($vep_version_data->{sub}) ? '.'.$vep_version_data->{sub} : ''),
      api_version   => $self->registry->software_version,
      input_headers => $self->get_Parser->headers,
    };

    if(my $mca = $self->get_adaptor('core', 'MetaContainer')) {
      $info->{db_name} = $mca->dbc->dbname;
      $info->{db_host} = $mca->dbc->host;
      $info->{db_version} = $mca->get_schema_version;
    }

    foreach my $as(@{$self->get_all_AnnotationSources}) {
      my $as_info = $as->info;
      $info->{version_data}->{$_} ||= $as_info->{$_} for grep {$_ ne 'custom_info'} keys %$as_info;
      $info->{cache_dir} ||= $as->dir if $as->can('dir');
      push @{$info->{custom_info}}, $as_info->{custom_info} if $as_info->{custom_info};
    }

    $self->{output_header_info} = $info;
  }

  return $self->{output_header_info};
}

sub valid_chromosomes {
  my $self = shift;

  if(!exists($self->{valid_chromosomes})) {

    my %valid = ();

    foreach my $as(@{$self->get_all_AnnotationSources}) {
      next unless $as->can('valid_chromosomes');
      $valid{$_} = 1 for @{$as->valid_chromosomes};
    }

    $self->{valid_chromosomes} = [sort keys %valid];
  }

  return $self->{valid_chromosomes};
}

# set some package variables to optimal values for speed
sub _set_package_variables {
  my $self = shift;

  # don't assert refs
  $self->{_assertions_bak} = $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
  $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0;

  # don't use rearrange
  $self->{_no_rearrange_bak} = $Bio::EnsEMBL::Utils::Argument::NO_REARRANGE;
  $Bio::EnsEMBL::Utils::Argument::NO_REARRANGE = 1;

  # avoid using transfer
  $self->{_no_transfer_bak} = $Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER;
  $Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER = 1;

  # suppress warnings that the FeatureAdpators spit if using no_slice_cache
  $self->{_verbose_bak} = Bio::EnsEMBL::Utils::Exception::verbose();
  Bio::EnsEMBL::Utils::Exception::verbose(1999);
}

sub _reset_package_variables {
  my $self = shift;

  $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = $self->{_assertions_bak};
  $Bio::EnsEMBL::Utils::Argument::NO_REARRANGE = $self->{_no_rearrange_bak};
  $Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER = $self->{_no_transfer_bak};

  Bio::EnsEMBL::Utils::Exception::verbose($self->{_verbose_bak});
}

1;