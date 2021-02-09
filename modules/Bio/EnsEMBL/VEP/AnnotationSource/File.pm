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
  'exact' => 1,
);


=head2 new

  Arg 1      : hashref $args
               {
                 config        => Bio::EnsEMBL::VEP::Config $config,
                 file          => string $filename,
                 format        => string $format (bed,bigwig,gff,gtf,vcf),
                 short_name    => (optional) string $short_name,
                 type          => (optional) string $type (overlap (default), exact),
                 report_coords => (optional) bool $report_coords,
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
  Description: Getter/setter for type of this source, either "overlap"
               or "exact".
                - overlap returns any source features that overlap input variants
                - exact requires that source features' coordinates match input variants exactly
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
      my ($vf_start, $vf_end) = ($vf->{start}, $vf->{end});
      if ($parser->seek($self->get_source_chr_name($chr), $vf_start - 1, $vf_end + 1)) {
        $parser->next();
      }

      # Different checks before annotating the VF:
      # - Check a record exist
      # - Check that the chromosomes names match between the input VF entry and the custom annotation file.
      #   There is a special case for the custom annotation file with chromosome names like 'chr1',
      #   as they are not present in the cache (chr_synonym.txt): we also check the 'raw' seqname with the source_chr_name.
      # - Check that the start of the custom annotation record is lower that the end of the input VF entry
      while($parser->{record} &&
            ($parser->get_seqname eq $self->get_source_chr_name($chr) || $parser->${get_raw_seqname} eq $self->get_source_chr_name($chr)) &&
            $parser->get_start <= $vf_end + 1) {
        $self->annotate_VariationFeature($vf);
        $parser->next();
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

  my $overlap_result = $self->_record_overlaps_VF($vf);

  push @{$vf->{_custom_annotations}->{$self->short_name}}, @{$self->_create_records($overlap_result)} if $overlap_result;
}


=head2 _create_records
 
  Example    : $records = $as->_create_records();
  Description: Create a custom annotation record from the current
               record as read from the annotation source.
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : annotate_VariationFeature()
  Status     : Stable

=cut

sub _create_records {
  return [{name => $_[0]->_get_record_name}];
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

  if($type eq 'overlap') {

    # account for insertions in Ensembl world where s = e+1
    my ($vs, $ve) = ($vf->{start}, $vf->{end});
    ($vs, $ve) = ($ve, $vs) if $vs > $ve;
    
    return overlap($parser->get_start, $parser->get_end, $vs, $ve);
  }
  elsif($type eq 'exact') {
    return $parser->get_start == $vf->{start} && $parser->get_end == $vf->{end};
  }
}

1;
