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

# EnsEMBL module for Bio::EnsEMBL::VEP::BaseRunner
#
#

=head1 NAME

Bio::EnsEMBL::VEP::BaseRunner - base runner class for VEP

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Base class shared by Bio::EnsEMBL::VEP::Runner and Bio::EnsEMBL::VEP::Haplo::Runner

=head1 METHODS

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


=head2 new

  Arg 1      : hashref $config
  Example    : $runner = Bio::EnsEMBL::VEP::Runner->new($config);
  Description: Creates a new runner object. The $config hash passed is
               used to create a Bio::EnsEMBL::VEP::Config object; see docs
               for this object and the vep script itself for allowed
               parameters.
  Returntype : Bio::EnsEMBL::VEP::BaseRunner
  Exceptions : throws on invalid configuration, see Bio::EnsEMBL::VEP::Config
  Caller     : vep
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # initialise self
  my $self = bless {}, $class;

  # get a config object
  $self->{_config} = Bio::EnsEMBL::VEP::Config->new(@_);

  return $self;
}


=head2 setup_db_connection

  Example    : $runner->setup_db_connection();
  Description: Sets up database connection. Also carries out a check of the
               assembly specified in the database versus one that the user
               may have supplied.
  Returntype : bool
  Exceptions : throws if:
                - assembly version from param does not match database
                - no assembly version found in database
  Caller     : vep
  Status     : Stable

=cut

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
      "\n You can use --offline to avoid a database connection entirely\n"
    ) if $config_assembly && $config_assembly ne $db_assembly;

    # update to database version
    $self->param('assembly', $db_assembly);

    if(!$self->param('assembly')) {
      throw("ERROR: No assembly version specified, use --assembly [version] or check the coord_system table in your core database\n");
    }
  }

  # update species, e.g. if user has input "human" we get "homo_sapiens"
  my $input_species = lc($self->param('species'));
  my $latin_species = $reg->get_alias($input_species);

  # fix for possible bug where a number is getting appended to the species name
  $latin_species =~ s/\d+$// if $latin_species =~ /^$input_species\d+$/;

  $self->species($latin_species);

  return 1;
}


=head2 get_all_AnnotationSources

  Example    : $sources = $runner->get_all_AnnotationSources();
  Description: Gets all AnnotationSources.
  Returntype : arrayref of Bio::EnsEMBL::VEP::AnnotationSource
  Exceptions : none
  Caller     : init(), _buffer_to_output()
  Status     : Stable

=cut

sub get_all_AnnotationSources {
  my $self = shift;

  if(!exists($self->{_annotation_sources})) {
    my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $self->config});
    $self->{_annotation_sources} = $asa->get_all;
  }

  return $self->{_annotation_sources};
}


=head2 get_output_file_handle

  Example    : $handle = $runner->get_output_file_handle();
  Description: Gets all file handle for writing output to.
  Returntype : glob
  Exceptions : throws if file exists or cannot write
  Caller     : run()
  Status     : Stable

=cut

sub get_output_file_handle {
  my $self = shift;

  unless(exists($self->{output_file_handle})) {

    my $output_file_handle = FileHandle->new();

    my $output_file_name = $self->param('output_file');
      
    # check if file exists
    if(-e $output_file_name && !$self->param('force_overwrite')) {
      throw("ERROR: Output file $output_file_name already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }
    
    # compress output?
    if(my $compress = $self->param('compress_output')) {

      # check required bin in path
      throw("ERROR: $compress not found in path\n") unless `which $compress 2>&1` =~ /\/$compress$/;

      if(uc($output_file_name) eq 'STDOUT') {
        $output_file_handle->open("| $compress -c |");
      }
      else {
        $output_file_handle->open("| $compress -c >$output_file_name") or throw("ERROR: Could not write to output file $output_file_name\n");
      }
    }
    # normal output
    else {
      if(uc($output_file_name) eq 'STDOUT') {
        $output_file_handle = *STDOUT;
      }
      else {
        $output_file_handle->open(">$output_file_name") or throw("ERROR: Could not write to output file $output_file_name\n");
      }
    }

    $self->{output_file_handle} = $output_file_handle;
  }

  return $self->{output_file_handle};
}


=head2 get_stats_file_handle

  Example    : $handle = $runner->get_stats_file_handle();
  Description: Gets all file handle for writing stats data to.
  Returntype : glob
  Exceptions : throws if file exists or cannot write
  Caller     : dump_stats()
  Status     : Stable

=cut

sub get_stats_file_handle {
  my $self = shift;
  my $type = shift;

  unless(exists($self->{stats_file_handle}->{$type})) {

    my $stats_file_root = $self->param('stats_file') || $self->param('output_file').'_summary';

    my $file_name;

    if($stats_file_root =~ m/^(.+)\.(.+?)$/) {
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
      throw("ERROR: Stats file $file_name already exists. Specify a different output file with --stats_file, overwrite existing file with --force_overwrite or disable stats with --no_stats\n");
    }
    
    my $fh = FileHandle->new();
    $fh->open(">$file_name") or throw("ERROR: Could not write to stats file $file_name\n");

    $self->{stats_file_handle}->{$type} = $fh;
  }

  return $self->{stats_file_handle}->{$type};
}


=head2 dump_stats

  Example    : $runner->dump_stats();
  Description: Writes all run stats to stats file
  Returntype : none
  Exceptions : none
  Caller     : finish()
  Status     : Stable

=cut

sub dump_stats {
  my $self = shift;

  unless($self->param('no_stats')) {

    if($self->param('stats_text')) {
      my $fh = $self->get_stats_file_handle('txt');
      $self->stats->dump_text($fh);
      close $fh;
    }

    if($self->param('stats_html')) {
      my $fh = $self->get_stats_file_handle('html');
      $self->stats->dump_html($fh);
      close $fh;
    }
  }
}


=head2 get_output_header_info

  Example    : $info = $runner->get_output_header_info();
  Description: Gets info hash including software version data, database info,
               input file headers, annotation source info
  Returntype : hashref
  Exceptions : none
  Caller     : init(), get_OutputFactory()
  Status     : Stable

=cut

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


=head2 valid_chromosomes

  Example    : $chrs = $runner->valid_chromosomes();
  Description: Compiles valid chromosomes across all annotation sources
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : get_Parser()
  Status     : Stable

=cut

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


=head2 _set_package_variables

  Example    : $runner->_set_package_variables();
  Description: Temporarily sets some package variables for speed
  Returntype : none
  Exceptions : none
  Caller     : next_output_line()
  Status     : Stable

=cut

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

  # HGVS shifting
  # Variable used when DB connection
  $self->{_shift_hgvs_db_bak} = $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME;
  # Variable used when offline
  $self->{_shift_hgvs_bak} = $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME;
  if (defined($self->param('shift_hgvs')) && $self->param('shift_hgvs') =~ /(0|1)/ ) {
    $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $self->param('shift_hgvs');
    $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $self->param('shift_hgvs');
  }

  # up/downstream distance
  if(my $distances = $self->param('distance')) {
    my ($u, $d) = @$distances;
    $d ||= $u;

    $self->{_upstream_bak} = $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE;
    $self->{_downstream_bak} = $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE;

    $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE = $u;
    $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE = $d;
  }
}


=head2 _reset_package_variables

  Example    : $runner->_reset_package_variables();
  Description: Re-sets package variables altered by _set_package_variables
  Returntype : none
  Exceptions : none
  Caller     : next_output_line()
  Status     : Stable

=cut

sub _reset_package_variables {
  my $self = shift;

  $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = $self->{_assertions_bak};
  $Bio::EnsEMBL::Utils::Argument::NO_REARRANGE = $self->{_no_rearrange_bak};
  $Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER = $self->{_no_transfer_bak};

  Bio::EnsEMBL::Utils::Exception::verbose($self->{_verbose_bak});

  $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $self->{_shift_hgvs_bak} if defined($self->{_shift_hgvs_bak});
  $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $self->{_shift_hgvs_db_bak} if defined($self->{_shift_hgvs_db_bak});

  $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE = $self->{_upstream_bak} if defined($self->{_upstream_bak});
  $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE = $self->{_downstream_bak} if defined($self->{_downstream_bak});
}

1;
