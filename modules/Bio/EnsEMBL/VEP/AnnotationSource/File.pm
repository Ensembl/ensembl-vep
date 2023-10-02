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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File - file-based annotation source

=head1 SYNOPSIS

# actually returns a Bio::EnsEMBL::VEP::AnnotationSource::File::VCF
my $as = Bio::EnsEMBL::VEP::AnnotationSource::File->new({
  config => $config,
  file   => "some_variants.vcf.gz",
  format => 'vcf',
  type   => 'exact'
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Base class for all custom file-based AnnotationSource classes.

Most child classes use ensembl-io parsers to read data from a custom annotation
file. This means they operate on a seek() and next() method of finding data in
the file.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);

our ($CAN_USE_TABIX_PM, $CAN_USE_BIGWIG);

BEGIN {
  if (eval q{ require Bio::DB::HTS::Tabix; 1 }) {
    $CAN_USE_TABIX_PM = 1;

    eval q{
      use Bio::EnsEMBL::VEP::AnnotationSource::File::BED;
      use Bio::EnsEMBL::VEP::AnnotationSource::File::VCF;
      use Bio::EnsEMBL::VEP::AnnotationSource::File::GFF;
      use Bio::EnsEMBL::VEP::AnnotationSource::File::GTF;
    };
  }

  if (eval q{ require Bio::DB::BigFile; 1 }) {
    $CAN_USE_BIGWIG = 1;

    eval q{
      use Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig;
    };
  }
}

my %FORMAT_MAP = (
  'vcf'     => 'VCF',
  'gff'     => 'GFF',
  'gtf'     => 'GTF',
  'bed'     => 'BED',
  'bigwig'  => 'BigWig',
);

my %VALID_TYPES = (
  'overlap' => 1,
  'within' => 1,
  'surrounding' => 1,
  'exact' => 1,
);


=head2 new

  Arg 1      : hashref $args
               {
                 config         => Bio::EnsEMBL::VEP::Config $config,
                 file           => string $filename,
                 format         => string $format (bed,bigwig,gff,gtf,vcf),
                 short_name     => (optional) string $short_name,
                 type           => (optional) string $type (overlap (default), within, surrounding, exact),
                 report_coords  => (optional) bool $report_coords,
                 overlap_cutoff => (optional) numeric $minimum_percentage_overlap (0 by default),
                 distance       => (optional) numeric $distance_to_overlapping_variant_ends (off by default),
                 same_type      => (optional) bool $only_match_identical_variant_classes (off by default),
                 reciprocal     => (optional) bool $calculate_reciprocal_overlap (off by default),
                 overlap_def    => (optional) string $overlap_definition (based on reciprocal by default),
                 num_records    => (optional) maximum number of records to show (50 by default),
                 summary_stats  => (optional) summary statistics: max, min, mean, count, sum
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::File->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::File object. Will
               actually return a sub-class depending on the specified format.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::File
  Exceptions : throws if:
                - no file given
                - invalid format specified
                - required module not installed (Bio::DB::HTS::Tabix or Bio::DB::BigFile)
  Caller     : AnnotationSourceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  throw("ERROR: No file given\n") unless $hashref->{file};
  $self->file($hashref->{file});

  $hashref->{short_name} = $self->short_name($hashref->{short_name} || (split '/', $self->file)[-1]);
  $hashref->{type} = $self->type($hashref->{type} || 'overlap');
  $self->report_coords(defined($hashref->{report_coords}) ? $hashref->{report_coords} : 0);

  $self->{overlap_cutoff} = $hashref->{overlap_cutoff} || 0;
  $self->{distance}       = $hashref->{distance};
  $self->{same_type}      = $hashref->{same_type}      || 0;
  $self->{reciprocal}     = $hashref->{reciprocal}     || 0;
  $self->{overlap_def}    = $hashref->{overlap_def};
  $self->{num_records}    = defined $hashref->{num_records} ? $hashref->{num_records} : 50;

  $self->{info} = { custom_info => $hashref };

  if(my $format = $hashref->{format}) {

    delete $hashref->{format};

    $format = lc($format);
    throw("ERROR: Unknown or unsupported format $format\n") unless $FORMAT_MAP{$format};

    if($format eq 'bigwig') {
      throw("ERROR: Cannot use format $format without Bio::DB::BigFile module installed\n") unless $CAN_USE_BIGWIG;
    }
    else {
      throw("ERROR: Cannot use format $format without Bio::DB::HTS::Tabix module installed\n") unless $CAN_USE_TABIX_PM;
    }
    $hashref->{_format} = $format;

    my $class = $self->module_prefix.'::AnnotationSource::File::'.$FORMAT_MAP{$format};
    eval "require $class";
    return $class->new({%$hashref, config => $self->config});
  }

  return $self;
}


=head2 file

  Arg 1      : (optional) string $file
  Example    : $file = $as->file()
  Description: Getter/setter for filename for this source.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub file {
  my $self = shift;
  $self->{file} = shift if @_;
  return $self->{file};
}


=head2 short_name

  Arg 1      : (optional) string $short_name
  Example    : $short_name = $as->short_name()
  Description: Getter/setter for short name for this source, used as the key
               in VEP output.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub short_name {
  my $self = shift;
  $self->{short_name} = shift if @_;
  return $self->{short_name};
}


=head2 type

  Arg 1      : (optional) string $type
  Example    : $type = $as->type()
  Description: Getter/setter for type of this source:
                - "overlap" returns source features that overlap input variants
                - "within" returns source features that are within the input variants
                - "surrounding" returns source features that completely surround the input variants
                - "exact" requires that source features' coordinates match input variants exactly
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;

  if(@_) {
    my $new_type = shift;
    throw("ERROR: New type \"$new_type\" is not valid\n") unless $VALID_TYPES{$new_type};
    $self->{type} = $new_type;
  }

  return $self->{type};
}


=head2 report_coords

  Arg 1      : (optional) bool $report_coords
  Example    : $report_coords = $as->report_coords()
  Description: Getter/setter for the report_coords parameter. If set to a true value,
               the coordinates of the source feature are returned in place of any
               identifier found in the source for the feature.
  Returntype : bool
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub report_coords {
  my $self = shift;
  $self->{report_coords} = shift if @_;
  return $self->{report_coords};
}


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $as->annotate_InputBuffer($ib);
  Description: Gets overlapping features from the source for the variants in
               the input buffer, and adds them as references to the  relevant
               VariationFeature.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my %by_chr;
  push @{$by_chr{$_->{chr}}}, $_ for @{$buffer->buffer};

  my $parser = $self->parser();

  # (Big)BED and (Big)Wig have the method 'get_raw_chrom' in ensembl-io
  # VCF and GXF (GFF and GTF) have the method 'get_raw_seqname' in ensembl-io
  my $get_raw_seqname = $parser->can('get_raw_seqname') ? 'get_raw_seqname' : 'get_raw_chrom';

  foreach my $chr(keys %by_chr) {
    foreach my $vf(@{$by_chr{$chr}}) {
      next if $vf->{vep_skip}; # avoid annotating previously skipped variants

      my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});
      if ($parser->seek($self->get_source_chr_name($chr), $vf_start - 1, $vf_end + 1)) {
        $parser->next();
      }

      my $stats = $self->{summary_stats};
      # Different checks before annotating the VF:
      # - Check a record exist
      # - Check that the chromosomes names match between the input VF entry and the custom annotation file.
      #   There is a special case for the custom annotation file with chromosome names like 'chr1',
      #   as they are not present in the cache (chr_synonym.txt): we also check the 'raw' seqname with the source_chr_name.
      # - Check that the start of the custom annotation record is lower that the end of the input VF entry
      my $record_count = 0;
      while($parser->{record} &&
            ($parser->get_seqname eq $self->get_source_chr_name($chr) || $parser->${get_raw_seqname} eq $self->get_source_chr_name($chr)) &&
            $parser->get_start <= $vf_end + 1) {
        # stop if exceeding the desired number of records and not calculating stats
        last if !defined $stats && $record_count > $self->{num_records};
        my $res = $self->annotate_VariationFeature($vf, $record_count);
        $record_count++ if $res;
        $parser->next();
      }

      # prepare statistics
      if (defined $stats) {
        my $annot_stats = $vf->{_custom_annotations_stats}->{$self->short_name};
        $annot_stats->{count} = $record_count if grep(/^count$/, @$stats);
        if ( grep(/^(mean)$/, @$stats) && defined $annot_stats->{sum} ) {
          $annot_stats->{mean} = $annot_stats->{sum} / $record_count;
        }
        delete $annot_stats->{sum} unless grep(/^sum$/, @$stats);
        $vf->{_custom_annotations_stats}->{$self->short_name} = $annot_stats;
      }
    }
  }
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
  return $self->{valid_chromosomes} ||= $self->parser->{tabix_file}->seqnames;
}


=head2 annotate_VariationFeature
 
  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature
  Arg 2      : $record_count
  Example    : $as->annotate_VariationFeature($vf);
  Description: Add custom annotations to the given variant using the
               current record as read from the annotation source.
  Returntype : none
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub annotate_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $record_count = shift;

  my ($overlap_result, $overlap_percentage) = $self->_record_overlaps_VF($vf);

  return unless ($overlap_result);

  my $stats      = $self->{summary_stats};
  my $get_scores = defined $stats;
  my $record = $self->_create_records($overlap_result, $get_scores);

  my $is_recorded = 0;
  if (@{$record}[0]->{'name'} =~  /^COSV/) {
    my ($matched_cosmic_record) = grep{$_->{'name'} eq @{$record}[0]->{'name'}} @{$vf->{_custom_annotations}->{$self->short_name}};
    if ($matched_cosmic_record){
      $is_recorded = 1;
      foreach my $key (keys %{@{$record}[0]->{'fields'}}) {
        unless (exists $matched_cosmic_record->{'fields'}->{$key}) {
          $matched_cosmic_record->{'fields'}{$key} = @{$record}[0]->{'fields'}->{$key};
        }   
      }
    }
  }

  if (!$is_recorded) {
    if (defined($self->{fields}) && grep {$_ eq "PC"} @{$self->{fields}}) {
      $record->[0]->{"fields"}->{"PC"} = $overlap_percentage;
    }

    # calculate summary statistics for custom annotation
    my $annot_stats = $vf->{_custom_annotations_stats}->{$self->short_name};
    my $value = $record->[0]->{score} if defined $stats;
    if (defined $value) {
      if ( grep(/^min$/, @$stats) ) {
        $annot_stats->{min} = $value if $value < ($annot_stats->{min} || '+inf');
      }
      if ( grep(/^max$/, @$stats) ) {
        $annot_stats->{max} = $value if $value > ($annot_stats->{max} || '-inf');
      }
      $annot_stats->{sum} += $value if grep(/^(sum|mean)$/, @$stats);
      $vf->{_custom_annotations_stats}->{$self->short_name} = $annot_stats;
    }

    if ( $record_count < $self->{num_records} ) {
      push @{$vf->{_custom_annotations}->{$self->short_name}}, @{$record};
    } elsif ( $record_count == $self->{num_records} ) {
      push @{$vf->{_custom_annotations}->{$self->short_name}}, { name => '...' };
    }
  }
  return 1;
}


=head2 _create_records
 
  Arg 1      : bool or hashref $overlap_result
  Arg 2      : bool $get_scores
  Example    : $records = $as->_create_records();
  Description: Create a custom annotation record from the current
               record as read from the annotation source.
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : annotate_VariationFeature()
  Status     : Stable

=cut

sub _create_records {
  my $self           = shift;
  my $overlap_result = shift;
  my $get_scores     = shift;

  my $record = [{ name  => $self->_get_record_name }];
  $record->[0]->{score} = $self->parser->get_score if $get_scores;
  return $record;
}


=head2 _get_record_name
 
  Example    : $record_name = $as->_get_record_name();
  Description: Get name for the current record using either ID as
               found in the source or record coordinates.
  Returntype : string
  Exceptions : none
  Caller     : _create_records()
  Status     : Stable

=cut

sub _get_record_name {
  my $self = shift;
  my $parser = $self->parser;

  my $name = $parser->get_name;

  return ($self->report_coords || !defined($name)) ?
    sprintf(
      '%s:%i-%i',
      $parser->get_seqname,
      $parser->get_start,
      $parser->get_end
    ) :
    $name;
}


=head2 _record_overlaps_VF
 
  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature
  Example    : $overlap_ok = $as->_record_overlaps_VF($vf);
  Description: Determine whether the given VariationFeature overlaps
               the current record, depending on the set type()
  Returntype : bool
  Exceptions : none
  Caller     : annotate_VariationFeature()
  Status     : Stable

=cut

sub _record_overlaps_VF {
  my $self = shift;
  my $vf = shift;

  my $parser = $self->parser();
  my $type = $self->type();
  my $overlap_cutoff = $self->{overlap_cutoff};
  my $distance = $self->{distance};
  my $reciprocal = $self->{reciprocal};

  my ($ref_start, $ref_end) = ($parser->get_start, $parser->get_end);
  $ref_start += 1 if defined $self->{_format} && $self->{_format} eq 'bigwig';

  if($type eq 'overlap' || $type eq 'within' || $type eq 'surrounding') {
    # account for insertions in Ensembl world where s = e+1
    my ($vs, $ve) = ($vf->{start}, $vf->{end});
    ($vs, $ve) = ($ve, $vs) if $vs > $ve;
    my $length = $ve - $vs + 1;

    # check if reference variant is within the input variant (if enabled)
    return 0 if $type eq "within" && ($vs > $ref_start || $ve < $ref_end);

    # check if reference variant completely surrounds the input variant (if enabled)
    return 0 if $type eq "surrounding" && ($vs < $ref_start || $ve > $ref_end);

    if (defined $distance) {
      return 0 unless abs($vs - $ref_start) <= $distance;
      return 0 unless abs($ve - $ref_end  ) <= $distance;
    }

    # check overlap percentage
    my @overlap_start = sort { $a <=> $b } ($vs, $ref_start);
    my @overlap_end   = sort { $a <=> $b } ($ve, $ref_end);
    my $overlap_percentage = 100 * (1 + $overlap_end[0] - $overlap_start[1]) / $length;

    return 0 if $overlap_percentage < $overlap_cutoff;

    if ($reciprocal) {
      # check bi-directional overlap - percentage of reference variant covered
      my $ref_length = $ref_end - $ref_start + 1;
      my $ref_overlap_percentage = 100 * (1 + $overlap_end[0] - $overlap_start[1]) / $ref_length;
      return 0 if $ref_overlap_percentage < $overlap_cutoff;

      # report minimum overlap
      if ($ref_overlap_percentage < $overlap_percentage) {
        $overlap_percentage = $ref_overlap_percentage;
      }
    }

    $overlap_percentage = sprintf("%.3f", $overlap_percentage);
    return overlap($ref_start, $ref_end, $vs, $ve), $overlap_percentage;
  }
  elsif($type eq 'exact') {
    my $match = $ref_start == $vf->{start} && $ref_end == $vf->{end};
    my $overlap_percentage = $match ? 100 : 0;
    return ( $match, $overlap_percentage );
  }
}

1;
