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

package Bio::EnsEMBL::VEP::Haplo::Runner;

use base qw(Bio::EnsEMBL::VEP::BaseRunner);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Haplo::TranscriptTree;
use Bio::EnsEMBL::VEP::Haplo::InputBuffer;
use Bio::EnsEMBL::VEP::Haplo::Parser::VCF;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->param('haplo', 1);

  return $self;
}

# dispatcher/runner for all initial setup from config
sub init {
  my $self = shift;

  return 1 if $self->{_initialized};

  # log start time
  $self->stats->start_time();

  # setup DB connection
  $self->setup_db_connection();

  # get chromosome synoyms
  $self->chromosome_synonyms($self->param('synonyms'));

  # get all annotation sources
  my $annotation_sources = $self->get_all_AnnotationSources();

  # setup FASTA file DB
  $self->fasta_db();

  my $buffer = $self->get_InputBuffer();

  $self->stats->info($self->get_output_header_info);

  return $self->{_initialized} = 1;
}

# run
sub run {
  my $self = shift;

  $self->init();

  my $input_buffer = $self->get_InputBuffer;

  my $count;
  
  while(my $vfs = $input_buffer->next()) {
    foreach my $as(@{$self->get_all_AnnotationSources}) {
      foreach my $thc(@{$as->annotate_InputBuffer($input_buffer)}) {
        my $tr = $thc->transcript;
        my $tr_stable_id = $tr->stable_id;

        return 1 if $count++ > 20;

        foreach my $ch(@{$thc->get_all_CDSHaplotypes}) {

          next if $ch->name eq $tr_stable_id.':REF';

          my $ph = $ch->get_ProteinHaplotype;
          my @out = (
            $tr_stable_id,
            $ch->name,
            $ph->name,
          );

          print join("\t", @out)."\n";
        }
      }
    }
  }

  # $self->dump_stats;

  return 1;
}

sub get_Parser {
  my $self = shift;

  if(!exists($self->{parser})) {

    # user given input data as string (REST)?
    if(my $input_data = $self->param('input_data')) {
      open IN, '<', \$input_data;
      $self->param('input_file', *IN);
    }

    $self->{parser} = Bio::EnsEMBL::VEP::Haplo::Parser::VCF->new({
      config            => $self->config,
      file              => $self->param('input_file'),
      valid_chromosomes => $self->get_valid_chromosomes,
    })
  }

  return $self->{parser};
}

sub get_InputBuffer {
  my $self = shift;

  if(!exists($self->{input_buffer})) {
    $self->{input_buffer} = Bio::EnsEMBL::VEP::Haplo::InputBuffer->new({
      config => $self->config,
      parser => $self->get_Parser,
      transcript_tree => $self->get_TranscriptTree
    });
  }

  return $self->{input_buffer};
}

sub get_TranscriptTree {
  my $self = shift;

  if(!exists($self->{transcript_tree})) {
    $self->{transcript_tree} = Bio::EnsEMBL::VEP::Haplo::TranscriptTree->new({
      config => $self->config,
      annotation_source => $self->get_all_AnnotationSources->[0]
    });
  }

  return $self->{transcript_tree};
}

1;