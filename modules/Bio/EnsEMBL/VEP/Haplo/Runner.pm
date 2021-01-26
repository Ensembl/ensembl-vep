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

# EnsEMBL module for Bio::EnsEMBL::VEP::Haplo::Runner
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::Runner - runner class for haplo

=head1 SYNOPSIS

my $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new($config_hashref);

$runner->run;

=head1 DESCRIPTION

The haplo runner class is used to run a Haplosaurus analysis given
the configuration provided to the new() method.

It should be the only class you need to instantiate manually if you
wish to incorporate Haplosaurus analysis into your code.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::Runner;

use base qw(Bio::EnsEMBL::VEP::BaseRunner);

use Scalar::Util qw(looks_like_number);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(get_compressed_filehandle);
use Bio::EnsEMBL::VEP::TranscriptTree;
use Bio::EnsEMBL::VEP::Haplo::InputBuffer;
use Bio::EnsEMBL::VEP::Haplo::Parser::VCF;

our $CAN_USE_JSON;

BEGIN {
  if(eval q{ use JSON; 1 }) {
    $CAN_USE_JSON = 1;
  }
}


=head2 new

  Arg 1      : hashref $config
  Example    : $runner = Bio::EnsEMBL::VEP::Haplo::Runner->new($config);
  Description: Creates a new haplo runner object. The $config hash passed is
               used to create a Bio::EnsEMBL::VEP::Config object; see docs
               for this object and the haplo script itself for allowed
               parameters.
  Returntype : Bio::EnsEMBL::VEP::Haplo::Runner
  Exceptions : throws on invalid configuration, see Bio::EnsEMBL::VEP::Config
  Caller     : haplo
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->param('haplo', 1);
  $self->param('output_file', 'haplo_output.txt') unless $self->config->_raw_config->{output_file};

  return $self;
}


=head2 init

  Example    : $runner->init()
  Description: Runs initialisation and setup for the runner. This includes:
                - setting up database connection
                - loading chromosome synonyms
                - getting AnnotationSource object
                - setting up FASTA file for sequence access
                - loading haplotype frequency data from disk
                - getting InputBuffer object
                - initialising stats
               Runs only once - subsequent calls are not executed.
  Returntype : bool
  Exceptions : none
  Caller     : run()
  Status     : Stable

=cut

sub init {
  my $self = shift;

  return 1 if $self->{_initialized};

  # check required modules etc
  $self->post_setup_checks();

  # log start time
  $self->stats->start_time();

  # setup DB connection
  $self->setup_db_connection();

  # get chromosome synoyms
  $self->chromosome_synonyms($self->param('synonyms'));

  # get annotation source
  $self->get_AnnotationSource();

  # setup FASTA file DB
  $self->fasta_db();

  # get haplotype frequency data
  $self->haplotype_frequencies($self->param('haplotype_frequencies'));

  my $buffer = $self->get_InputBuffer();

  $self->stats->info($self->get_output_header_info);

  return $self->{_initialized} = 1;
}


=head2 post_setup_checks

  Example    : $runner->post_setup_checks();
  Description: Checks things from configuration.
  Returntype : none
  Exceptions : none
  Caller     : init()
  Status     : Stable

=cut


sub post_setup_checks {
  my $self = shift;

  die("ERROR: JSON module not installed\n") if $self->param('json') && !$CAN_USE_JSON;
}


=head2 run

  Example    : $runner->run();
  Description: Runs Haplosaurus using the config as set up with new(). Writes
               output to file/handle.
  Returntype : bool
  Exceptions : none
  Caller     : haplo
  Status     : Stable

=cut

sub run {
  my $self = shift;

  $self->init();

  my $input_buffer = $self->get_InputBuffer;
  my $as = $self->get_AnnotationSource;

  $self->_set_package_variables();

  my $vfs = $input_buffer->next();
  
  while(@$vfs) {
    $self->dump_TranscriptHaplotypeContainer($_) for @{$as->annotate_InputBuffer($input_buffer)};
    $vfs = $input_buffer->next();
  }

  $self->_reset_package_variables();

  # Send a warning message if the output is empty (i.e. no TranscriptHaplotypeContainer)
  if (!$self->{_output_lines_count} || $self->{_output_lines_count} == 0) {
    warning("Haplosaurus can't find transcripts overlapping your variant(s). The output is empty.");
  }

  return 1;
}


=head2 dump_TranscriptHaplotypeContainer

  Arg 1      : Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer $thc
  Example    : $runner->dump_TranscriptHaplotypeContainer($thc);
  Description: Writes a TranscriptHaplotypeContainer object to a flat file
               ($self->get_output_file_handle), either in tab-delimited or
               JSON format.
  Returntype : none
  Exceptions : none
  Caller     : run()
  Status     : Stable

=cut

sub dump_TranscriptHaplotypeContainer {
  my ($self, $thc) = @_;

  my $fh = $self->get_output_file_handle;

  if($self->param('json')) {
    my $json = $self->{_json} || JSON->new->convert_blessed;

    $thc->_dont_export($_) for @{$self->param('dont_export') || []};

    print $fh $json->encode($thc);
    print $fh "\n";

    $self->{_output_lines_count} ++;

    return;
  }

  my $tr = $thc->transcript;
  my $tr_stable_id = $tr->stable_id;

  foreach my $ch(
    sort {$b->count <=> $a->count || $a->name cmp $b->name}
    grep {!$_->is_reference}
    @{$thc->get_all_CDSHaplotypes}
  ) {

    my $sample_counts = $ch->get_all_sample_counts;

    my $freqs = "";
    my $ph = $ch->get_ProteinHaplotype;

    if(my $freq_data = $self->haplotype_frequencies->{$ph->_hex}) {
      my $tr_freq_data;

      if(scalar @$freq_data > 1) {
        if($freq_data->[0]->{transcript}) {
          ($tr_freq_data) = grep {$_->{transcript} eq $tr_stable_id} @$freq_data;
        }
        else {
          foreach my $tr_freq_data(@$freq_data) {
            $freqs = join(",", map {$_.'='.$tr_freq_data->{$_}} keys %$tr_freq_data);
          }
        }
      }
      else {
        $tr_freq_data = $freq_data->[0];
      }

      $freqs = join(",",
        map {$_.'='.$tr_freq_data->{$_}}
        grep {$_ ne 'transcript'}
        keys %$tr_freq_data
      );
    }

    my @out = (
      $tr_stable_id,
      $ch->name,
      join(",", @{$ch->get_all_flags}),
      $ph->name,
      join(",", @{$ph->get_all_flags}),
      $freqs,
      join(",", map {$_->variation_name} @{$ch->get_all_VariationFeatures}),
      join(",", map {$_.':'.$sample_counts->{$_}} sort keys %$sample_counts)
    );

    print $fh join("\t", @out)."\n";

    $self->{_output_lines_count} ++;
  }
}


=head2 get_AnnotationSource

  Example    : my $as = $runner->get_AnnotationSource();
  Description: Gets the AnnotationSource configured for this runner.
  Returntype : Bio::EnsEMBL::VEP::Haplo::AnnotationType::Transcript
  Exceptions : none
  Caller     : init(), run(), get_TranscriptTree()
  Status     : Stable

=cut

sub get_AnnotationSource {
  return $_[0]->get_all_AnnotationSources->[0];
}


=head2 get_Parser

  Example    : my $parser = $runner->get_Parser();
  Description: Gets the Parser configured for this runner. Currently Haplosaurus
               supports only VCF input.
  Returntype : Bio::EnsEMBL::VEP::Haplo::Parser::VCF
  Exceptions : none
  Caller     : getInputBuffer()
  Status     : Stable

=cut

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
      valid_chromosomes => $self->valid_chromosomes,
    })
  }

  return $self->{parser};
}


=head2 get_InputBuffer

  Example    : my $ib = $runner->get_InputBuffer();
  Description: Gets the InputBuffer configured for this runner.
  Returntype : Bio::EnsEMBL::VEP::Haplo::InputBuffer
  Exceptions : none
  Caller     : init(), run()
  Status     : Stable

=cut

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


=head2 get_TranscriptTree

  Example    : my $parser = $runner->get_TranscriptTree();
  Description: Gets the TranscriptTree configured for this runner.
  Returntype : Bio::EnsEMBL::VEP::TranscriptTree
  Exceptions : none
  Caller     : get_InputBuffer()
  Status     : Stable

=cut

sub get_TranscriptTree {
  my $self = shift;

  if(!exists($self->{transcript_tree})) {
    $self->{transcript_tree} = Bio::EnsEMBL::VEP::TranscriptTree->new({
      config => $self->config,
      annotation_source => $self->get_AnnotationSource
    });
  }

  return $self->{transcript_tree};
}


=head2 haplotype_frequencies
  
  Arg 1      : string $filename
  Example    : my $freqs = $runner->haplotype_frequencies($freqs_file);
  Description: Gets a hashref of haplotype frequencies keyed on the md5 hex
               of the haplotype's sequence. If a filename is provided, the
               frequencies are loaded from this file. The file should contain
               a header line starting with "#" indicating the names of the
               columns; all columns except those named "hex" or "transcript"
               are treated as representing a population/grouping name.
  Returntype : hashref
  Exceptions : throws if no hex field in header or invalid data in file
  Caller     : init(), dump_TranscriptHaplotypeContainer()
  Status     : Stable

=cut

sub haplotype_frequencies {
  my ($self, $file) = @_;

  $self->{_haplotype_frequencies} ||= {};

  if($file) {
    my $fh;

    if(-B $file) {
      $fh = get_compressed_filehandle($file);
    }
    else {
      open $fh, $file or throw("ERROR: $!");
    }

    my @headers;

    while(<$fh>) {
      chomp;

      if(/^\#/) {
        s/^\#+\s*//;
        @headers = split;

        throw("ERROR: No hex header found in haplotype frequency file $file\n") unless grep {$_ eq 'hex'} @headers;
      }
      else {
        my @data = split;

        # no headers? we can have a look at the file to see if we can guess the format
        unless(@headers) {
          
          # first col should be hex, if not then we can't do much else
          if($data[0] =~ /^[a-fA-F0-9]{32}$/) {
            $headers[0] = 'hex';
          }
          else {
            throw("ERROR: Cannot parse haplotype frequency file $file, first column should be an md5sum\n");
          }

          # examine remaining columns to see if they look like numbers
          my $f = 1;
          my $u = 1;

          for my $i(1..$#data) {
            $headers[$i] = (looks_like_number($data[$i]) && $data[$i] >= 0 && $data[$i] <= 1) ? 'freq'.$f++ : 'unknown'.$u++;
          }

          throw("ERROR: Could not find any data that looks like frequencies in haplotype frequency file $file\n") unless $f > 1;
        }

        throw("ERROR: Number of columns (".(scalar @data).") does not match headers (".(scalar @headers).")\n") unless scalar @data == scalar @headers;

        my %row = map {$headers[$_] => $data[$_]} 0..$#data;
        my $hex = delete $row{hex};

        push @{$self->{_haplotype_frequencies}->{$hex}}, \%row;
      }
    }

    close $fh;
  }

  return $self->{_haplotype_frequencies};
}

1;
