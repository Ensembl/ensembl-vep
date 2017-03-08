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

use Scalar::Util qw(looks_like_number);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(get_compressed_filehandle);
use Bio::EnsEMBL::VEP::TranscriptTree;
use Bio::EnsEMBL::VEP::Haplo::InputBuffer;
use Bio::EnsEMBL::VEP::Haplo::Parser::VCF;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->param('haplo', 1);
  $self->param('output_file', 'haplo_output.txt') unless $self->config->_raw_config->{output_file};

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

# run
sub run {
  my $self = shift;

  $self->init();

  my $input_buffer = $self->get_InputBuffer;
  my $as = $self->get_AnnotationSource;

  my $count;

  $self->_set_package_variables();

  my $vfs = $input_buffer->next();
  
  while(@$vfs) {
    $self->dump_TranscriptHaplotypeContainer($_) for @{$as->annotate_InputBuffer($input_buffer)};
    $vfs = $input_buffer->next();
  }

  $self->_reset_package_variables();

  # $self->dump_stats;

  return 1;
}

sub dump_TranscriptHaplotypeContainer {
  my ($self, $thc) = @_;

  my $fh = $self->get_output_file_handle;

  my $tr = $thc->transcript;
  my $tr_stable_id = $tr->stable_id;

  foreach my $ch(
    sort {$b->count <=> $a->count || $a->name cmp $b->name}
    grep {!$_->is_reference}
    @{$thc->get_all_CDSHaplotypes}
  ) {

    my $sample_counts = $ch->get_all_sample_counts;

    foreach my $sample(keys %$sample_counts) {

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
        $sample,
        $sample_counts->{$sample},
      );

      print $fh join("\t", @out)."\n";
    }
  }
}

sub get_AnnotationSource {
  return $_[0]->get_all_AnnotationSources->[0];
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
      valid_chromosomes => $self->valid_chromosomes,
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
    $self->{transcript_tree} = Bio::EnsEMBL::VEP::TranscriptTree->new({
      config => $self->config,
      annotation_source => $self->get_all_AnnotationSources->[0]
    });
  }

  return $self->{transcript_tree};
}

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