=head1 LICENSE

Copyright [2016-2026] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig - BigWig annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig->new({
  config => $config,
  file   => "my_scores.bw",
  type   => "exact"
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

BigWig format custom file annotation source.

BigWig format is used for storing dense continuous data e.g. per-base scores.

https://genome.ucsc.edu/goldenpath/help/bigWig.html

This module requires Bio::DB::BigFile to be installed.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::IO::Parser::BigWig;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File);

# Only min/max are served from the zoom reductions - a bin's minVal/maxVal are
# exact extrema over the bin. sum/count/mean can't be read back from zoom (they
# come out inflated: kent bbiRead.c normalizeCount rounds each bin's apportioned
# base count up to an integer), and they force an O(span) walk anyway, so
# annotate_InputBuffer sends any request touching them to the base class.
#
# Window sizes for the decomposition, as methods so a subclass (or test) can
# override them:
#   zoom_grid     interior zoom bin width (G)
#   exact_end     exact per-base window at each end (E); the E > G invariant caps
#                 any zoom bin straddling the interior boundary inside an exact end
#   zoom_min_span spans <= 2E skip zoom and scan per-base (unchanged small-variant path)
sub zoom_grid     { 16384 }
sub exact_end     { my $self = shift; 2 * $self->zoom_grid }
sub zoom_min_span { my $self = shift; 2 * $self->exact_end }


=head2 parser

  Example    : $parser = $as->parser();
  Description: Get ensembl-io parser to read from file
  Returntype : Bio::EnsEMBL::IO::Parser::BigWig. Note this includes
               a hacky implementation using fork() to avoid issues with the kent
               source tree calling exit() if a file doesn't exist,
               causing the whole script to exit without explanation.
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;

  unless(exists($self->{parser})) {

    # This is something of a hack.
    # The BigWig parser uses the Kent source tree which calls exit() if the file doesn't exist.
    # We can't use -e $self->file as the file might be a remote one.
    # So, we fork, attempt to open the file, and check the return code of the forked process.
    my $pid = fork();

    # parent process, capture the return code and die nicely if it failed
    if($pid) {
      if(waitpid($pid, 0) > 0) {
        my $rc = $? >> 8;
        throw("ERROR: Failed to open ".$self->file."\n") if $rc > 0;
      }
    }

    # child process, attempt to open the file
    else {
      my $test = Bio::EnsEMBL::IO::Parser::BigWig->open($self->file);
      exit(0);
    }

    $self->{parser} = Bio::EnsEMBL::IO::Parser::BigWig->open($self->file);
  }
  return $self->{parser};
}


=head2 valid_chromosomes

  Example    : $chrs = $as->valid_chromosomes();
  Description: Gets valid chromosome names for this source
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub valid_chromosomes {
  my $self = shift;
  return $self->{valid_chromosomes} ||= [keys %{$self->parser->{cache}->{chromosomes}}];
}


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $as->annotate_InputBuffer($ib);
  Description: BigWig-specific override of the base annotation loop. When
               summary_stats requests only min and/or max, computing them over a
               large reference span (e.g. a multi-Mb structural variant) with the
               default per-base walk is O(span) in time and memory and can
               effectively never finish. For such spans this override reads the
               interior from the bigWig's precomputed zoom reductions instead
               (see _span_min_max), which is orders of magnitude faster and keeps
               min/max bit-exact. Small variants (SNVs, indels, small SVs), any
               summary_stats combination that also needs sum/mean/count, and every
               non-summary case fall through to the base class, so their behaviour
               is unchanged.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my $stats = $self->{summary_stats};

  return $self->SUPER::annotate_InputBuffer($buffer)
    unless defined($stats) && $self->_can_use_zoom_summary($stats);

  my @exact_vfs;
  foreach my $vf (@{$buffer->buffer}) {
    if ($vf->{vep_skip} || $self->_span_length($vf) <= $self->zoom_min_span) {
      push @exact_vfs, $vf;
    }
    else {
      $self->_annotate_large_span_vf($vf, $stats);
    }
  }

  # run the unchanged base-class loop for the rest by pointing it at just them
  if (@exact_vfs) {
    my $all = $buffer->buffer;
    my $ok = eval {
      $buffer->buffer(\@exact_vfs);
      $self->SUPER::annotate_InputBuffer($buffer);
      1;
    };
    $self->SUPER::annotate_InputBuffer($buffer);
    $buffer->buffer($all);
  }

  return;
}


=head2 _span_length

  Description: Reference-span length (bp) of a VariationFeature. Insertions
               (start = end + 1 in Ensembl coordinates) give 0.
  Returntype : integer

=cut

sub _span_length {
  my ($self, $vf) = @_;
  my ($s, $e) = ($vf->{start}, $vf->{end});
  ($s, $e) = ($e, $s) if $s > $e;
  return $e - $s + 1;
}


=head2 _can_use_zoom_summary

  Arg 1      : arrayref $summary_stats
  Description: True only when the requested statistics are a non-empty subset of
               min/max and the configured matching options are the plain defaults
               that the zoom decomposition can reproduce (type "overlap", no
               overlap cutoff / distance / reciprocal / same_type restriction,
               and no per-record PC field). Anything else falls back to the
               exact per-base path.
  Returntype : bool

=cut

sub _can_use_zoom_summary {
  my ($self, $stats) = @_;
  return 0 unless @$stats && !grep { $_ ne 'min' && $_ ne 'max' } @$stats;
  return 0 unless $self->type eq 'overlap';
  return 0 if $self->{overlap_cutoff};
  return 0 if defined $self->{distance};
  return 0 if $self->{reciprocal};
  return 0 if $self->{same_type};
  return 0 if $self->{fields} && grep { $_ eq 'PC' } @{$self->{fields}};
  return 1;
}


=head2 _annotate_large_span_vf

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature
  Arg 2      : arrayref $summary_stats  (a subset of min/max)
  Description: Compute and attach the requested min/max statistics for a single
               large-span variant using the aligned exact/zoom decomposition, and
               populate the same truncated per-record list the base class would
               produce so the output shape is unchanged.
  Returntype : none

=cut

sub _annotate_large_span_vf {
  my ($self, $vf, $stats) = @_;

  my $parser = $self->parser;
  my $chr    = $self->get_source_chr_name($vf->{chr});
  my $seq_id = $parser->cache->{chromosomes}->{$chr};
  return unless defined $seq_id;    # chromosome not present in this bigWig

  # convert Ensembl 1-based inclusive [start, end] to 0-based half-open [s, e),
  # clamped to the contig - matching the set of intervals the base class scans
  my ($vs, $ve) = ($vf->{start}, $vf->{end});
  ($vs, $ve) = ($ve, $vs) if $vs > $ve;
  my $s = $vs - 1;
  my $e = $ve;
  $s = 0 if $s < 0;
  my $chr_len = $parser->cache->{chr_sizes}->{$chr};
  $e = $chr_len if defined($chr_len) && $e > $chr_len;
  return if $e <= $s;

  my $fh = $parser->open_file;
  my ($min, $max) = $self->_span_min_max($fh, $seq_id, $s, $e);
  return unless defined $max;    # undef only when [s, e) has no data

  my $annot = $vf->{_custom_annotations_stats}->{$self->short_name} ||= {};
  $annot->{min} = $min if grep { $_ eq 'min' } @$stats;
  $annot->{max} = $max if grep { $_ eq 'max' } @$stats;

  # populate the per-record list (up to num_records, then '...') for output
  # parity with the base class; stats are only emitted when this key exists
  $self->_collect_span_records($vf, $fh, $seq_id, $s, $e);

  return;
}


=head2 _span_min_max

  Arg 1      : Bio::DB::bbiFile $fh
  Arg 2      : string $seq_id  (bigWig-internal chromosome name)
  Arg 3      : integer $s      (0-based half-open start)
  Arg 4      : integer $e      (0-based half-open end)
  Description: (min, max) over [s, e), bit-exact vs the per-base walk: the two ends
               (width exact_end, and any span up to 2*exact_end) are read exactly
               with bigWigIntervalQuery, the interior extrema come from the zoom
               reductions via bigWigSummaryArrayExtended on a grid of width
               zoom_grid.
  Returntype : list (min, max); both undef if [s, e) has no data

=cut

sub _span_min_max {
  my ($self, $fh, $seq_id, $s, $e) = @_;

  my ($min, $max);
  my $see = sub {
    my ($lo_val, $hi_val) = @_;
    $min = $lo_val if !defined($min) || $lo_val < $min;
    $max = $hi_val if !defined($max) || $hi_val > $max;
  };

  my $exact = sub {
    my ($ws, $we) = @_;
    return if $we <= $ws;
    my $list = $fh->bigWigIntervalQuery("$seq_id", $ws, $we);
    for (my $i = $list->head; $i; $i = $i->next) {
      $see->($i->value, $i->value);
    }
  };

  my $E = $self->exact_end;
  if ($e - $s <= 2 * $E) {
    $exact->($s, $e);
    return ($min, $max);
  }

  my ($lo, $hi) = ($s + $E, $e - $E);
  $exact->($s, $lo);
  $exact->($hi, $e);

  my $nbins = int(($hi - $lo) / $self->zoom_grid);
  $nbins = 1 if $nbins < 1;
  my $ext = $fh->bigWigSummaryArrayExtended("$seq_id", $lo, $hi, $nbins);
  if (defined $ext) {
    foreach my $bin (@$ext) {
      next unless $bin->{validCount} > 0;
      $see->($bin->{minVal}, $bin->{maxVal});
    }
  }

  return ($min, $max);
}


=head2 _collect_span_records

  Description: Fill $vf->{_custom_annotations} with up to num_records records
               (then a '...' marker) taken from the start of the span, matching
               the truncated list the base class emits. Only the left end of the
               span is scanned so this stays cheap on huge spans; for the dense
               per-base tracks this fast path targets these are exactly the first
               records the base class would report.
  Returntype : none

=cut

sub _collect_span_records {
  my ($self, $vf, $fh, $seq_id, $s, $e) = @_;

  my $num_records = $self->{num_records};

  # ensure the key exists (even with num_records 0) so the stats are emitted
  my $records = $vf->{_custom_annotations}->{$self->short_name} ||= [];

  if ($num_records <= 0) {
    push @$records, { name => '...' };
    return;
  }

  my $we = $s + $self->exact_end;
  $we = $e if $e < $we;
  my $list = $fh->bigWigIntervalQuery("$seq_id", $s, $we);

  for (my $i = $list->head; $i; $i = $i->next) {
    if (scalar(@$records) >= $num_records) {
      push @$records, { name => '...' };
      return;
    }
    my $name = $self->report_coords
      ? sprintf('%s:%i-%i', $seq_id, $i->start + 1, $i->end)
      : $i->value;
    push @$records, { name => $name };
  }

  # the span continues past the scanned window, so mirror the base class marker
  push @$records, { name => '...' } if $e > $we;

  return;
}


=head2 _get_record_name

  Example    : $record_name = $as->_get_record_name();
  Description: Get name for the current record using either score as
               found in the bigWig or record coordinates.
  Returntype : string
  Exceptions : none
  Caller     : _create_records()
  Status     : Stable

=cut

sub _get_record_name {
  my $self = shift;
  my $parser = $self->parser;

  return $self->report_coords ?
    sprintf(
      '%s:%i-%i',
      $parser->get_seqname,
      $parser->get_start,
      $parser->get_end
    ) :
    $parser->get_score;
}

1;