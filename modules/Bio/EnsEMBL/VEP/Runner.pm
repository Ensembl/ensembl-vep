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

# EnsEMBL module for Bio::EnsEMBL::VEP::Runner
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Runner - runner class for VEP

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Runner;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::VEP::Utils qw(get_time);
use Bio::EnsEMBL::VEP::Constants;
use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::Parser;
use Bio::EnsEMBL::VEP::InputBuffer;
use Bio::EnsEMBL::VEP::OutputFactory;
use Bio::EnsEMBL::VEP::AnnotationSourceAdaptor;

# don't assert refs
$Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0;

# don't use rearrange
$Bio::EnsEMBL::Utils::Argument::NO_REARRANGE = 1;

# avoid using transfer
$Bio::EnsEMBL::Variation::TranscriptVariationAllele::NO_TRANSFER = 1;

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

# dispatcher/runner for all initial setup from config
sub init {
  my $self = shift;

  return 1 if $self->{_initialized};

  # setup DB connection
  $self->setup_db_connection();

  # get all annotation sources
  my $annotation_sources = $self->get_all_AnnotationSources();

  # setup FASTA file DB
  my $fasta_db = $self->setup_fasta_db();

  my $buffer = $self->get_InputBuffer();

  return $self->{_initialized} = 1;
}

# run
# sub run {
#   my $self = shift;

#   $self->init();

#   my $input_buffer = $self->input_buffer;

#   while(my $vfs = $input_buffer->next()) {
#     last unless scalar @$vfs;

#     foreach my $as(@{$self->get_all_AnnotationSources}) {
#       $as->annotate_InputBuffer($input_buffer);
#     }
#   }
# }

sub next_output_line {
  my $self = shift;

  my $output_buffer = $self->{_output_buffer} ||= [];

  return shift @$output_buffer if @$output_buffer;

  $self->init();

  my $input_buffer = $self->get_InputBuffer;
  my $vfs = $input_buffer->next();

  if($vfs && scalar @$vfs) {
    my $output_factory = $self->get_OutputFactory;

    foreach my $as(@{$self->get_all_AnnotationSources}) {
      $as->annotate_InputBuffer($input_buffer);
    }
      
    $input_buffer->finish_annotation;

    push @$output_buffer, @{$output_factory->get_all_lines_by_InputBuffer($input_buffer)};
  }

  return @$output_buffer ? shift @$output_buffer : undef;
}

sub setup_db_connection {
  my $self = shift;

  return if $self->param('offline');

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
  $self->species($reg->get_alias($self->param('species')));

  return 1;
}

sub setup_fasta_db {
  my $self = shift;

  if(!exists($self->config->{_fasta_db})) {
    my $fasta_db;

    if(my $fasta_file = $self->param('fasta')) {

      $fasta_db = setup_fasta(
        -FASTA => $fasta_file,
        -ASSEMBLY => $self->param('assembly'),
        -OFFLINE => $self->param('offline'),
      );
    }

    $self->config->{_fasta_db} = $fasta_db;
  }

  return $self->config->{_fasta_db};
}

sub get_all_AnnotationSources {
  my $self = shift;

  if(!exists($self->{_annotation_sources})) {
    my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $self->config});
    $self->{_annotation_sources} = $asa->get_all;
  }

  return $self->{_annotation_sources};
}

sub get_Parser {
  my $self = shift;

  if(!exists($self->{parser})) {

    # user given input data as string (REST)?
    if(my $input_data = $self->param('input_data')) {
      open IN, '<', \$input_data;
      $self->param('input_file', *IN);
    }

    $self->{parser} = Bio::EnsEMBL::VEP::Parser->new({
      config => $self->config,
      format => $self->param('format'),
      file   => $self->param('input_file')
    })
  }

  return $self->{parser};
}

sub get_InputBuffer {
  my $self = shift;

  if(!exists($self->{input_buffer})) {
    $self->{input_buffer} = Bio::EnsEMBL::VEP::InputBuffer->new({
      config => $self->config,
      parser => $self->get_Parser
    });
  }

  return $self->{input_buffer};
}

sub get_OutputFactory {
  my $self = shift;

  if(!exists($self->{output_factory})) {
    $self->{output_factory} = Bio::EnsEMBL::VEP::OutputFactory->new({
      config      => $self->config,
      format      => $self->param('output_format'),
      header_info => $self->get_output_header_info,
    });
  }

  return $self->{output_factory};
}

sub get_output_header_info {
  my $self = shift;

  if(!exists($self->{output_header_info})) {

    my $info = {
      time          => get_time,
      vep_version   => $Bio::EnsEMBL::VEP::Constants::VERSION,
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
      $info->{version_data}->{$_} ||= $as_info->{$_} for keys %$as_info;
      $info->{cache_dir} ||= $as->dir if $as->can('dir');
    }

    $self->{output_header_info} = $info;
  }

  return $self->{output_header_info};
}

1;