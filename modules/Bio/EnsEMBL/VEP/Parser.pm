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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser - base class for input parsers

=head1 SYNOPSIS

# note this actually returns a Bio::EnsEMBL::VEP::Parser::VCF
my $parser = Bio::EnsEMBL::VEP::Parser->new({
  config => $config,
  file   => 'variants.vcf',
  format => 'vcf'
});

my $vf = $parser->next();

=head1 DESCRIPTION

Base class for Parsers; typically a Parser sub-class will
read data from an ensembl-io parser class.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Parser;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(get_compressed_filehandle);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(trim_sequences);

use Bio::EnsEMBL::VEP::Parser::VCF;
use Bio::EnsEMBL::VEP::Parser::VEP_input;
use Bio::EnsEMBL::VEP::Parser::ID;
use Bio::EnsEMBL::VEP::Parser::HGVS;
use Bio::EnsEMBL::VEP::Parser::Region;
use Bio::EnsEMBL::VEP::Parser::SPDI;

use Scalar::Util qw(openhandle looks_like_number);
use FileHandle;

my %FORMAT_MAP = (
  'vcf'     => 'VCF',
  'ensembl' => 'VEP_input',
  'id'      => 'ID',
  'hgvs'    => 'HGVS',
  'region'  => 'Region',
  'spdi'    => 'SPDI'
);


=head2 new

  Arg 1      : hashref $args
               {
                 config    => Bio::EnsEMBL::VEP::Config,
                 file      => string or filehandle,
                 format    => (optional) string (vcf, ensembl, id, hgvs),
                 delimiter => (optional) string,
               }
  Example    : $parser = Bio::EnsEMBL::VEP::Parser->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::Parser object. Will
               actually return a sub-class depending on either the specified
               or detected format.
  Returntype : Bio::EnsEMBL::VEP::Parser
  Exceptions : throws if:
                - no file given
                - invalid format specified
                - unable to detect format
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  $self->add_shortcuts([qw(dont_skip check_ref chr lrg minimal delimiter lookup_ref)]);

  my $hashref = $_[0];

  die("ERROR: No input file given\n") unless $hashref->{file};
  $self->file($hashref->{file});

  $self->delimiter($hashref->{delimiter}) if $hashref->{delimiter};

  $self->line_number(0);

  $self->valid_chromosomes({map {$_ => 1} @{$hashref->{valid_chromosomes} || []}});

  if(my $format = $hashref->{format}) {

    delete $hashref->{format};

    # detect format
    if(lc($format) eq 'guess' || lc($format) eq 'detect' || lc($format) eq 'auto') {
      $format = $self->detect_format();
      $self->warning_msg("No input file format specified - detected $format format") if $self->param('verbose') && defined $format;
    }

    die("ERROR: Can't detect input format\n") unless $format;

    $format = lc($format);
    die("ERROR: Unknown or unsupported input format '$format'\n") unless $FORMAT_MAP{$format};

    $self->param('format', $format);

    my $class = 'Bio::EnsEMBL::VEP::Parser::'.$FORMAT_MAP{$format};
    return $class->new({%$hashref, config => $self->config, delimiter => $self->delimiter});
  }

  return $self;
}


=head2 next

  Example    : $vf = $parser->next();
  Description: Fetches the next valid VariationFeature from the input file
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : none
  Caller     : InputBuffer
  Status     : Stable

=cut

sub next {
  my $self = shift;

  my $cache = $self->{_vf_cache} ||= [];

  if(!scalar @$cache) {
    push @$cache, @{$self->create_VariationFeatures()};
  }

  my $vf = shift @$cache;
  return $vf unless $vf;

  unless($self->validate_vf($vf) || $self->{dont_skip}) {
    return $self->next();
  }

  return $vf;
}


=head2 headers

  Example    : $headers = $parser->headers();
  Description: Gets headers from the input file. Just a stub here.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub headers {
  return [];
}


=head2 file

  Arg 1      : (optional) string $filename or glob *filehandle
  Example    : $file = $parser->file();
  Description: Getter/setter for input file. This may be an open filehandle,
               a file name (optionally gzip or bgzip compressed), or the string
               'STDIN' to represent reading from STDIN.
  Returntype : string or filehandle
  Exceptions : throws if filename does not exist
  Caller     : new(), parser()
  Status     : Stable

=cut

sub file {
  my $self = shift;

  if(@_) {
    my $file = shift;

    $self->{file_bak} = $file;

    if(uc($file) eq 'STDIN') {
      $file = *STDIN;
    }
    elsif(!openhandle($file)) {
      throw("ERROR: File \"$file\" does not exist\n") unless -e $file;

      # file is compressed?
      $file = get_compressed_filehandle($file, 1) if -f $file && -B $file;
    }

    $self->{file} = $file;
  }

  return $self->{file};
}


=head2 skip_empty_lines

  Example    : $num_skipped = $parser->skip_empty_lines();
  Description: Skip over any empty lines that may be in the input,
               returning the number skipped.
               Not used by all parsers, so not executed in the next()
               method in this base class.
  Returntype : int
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub skip_empty_lines {
  my $self = shift;
  my $parser = $self->parser;

  my $skipped = 0;

  # allow for empty lines
  while(defined($parser->{record}) && $parser->{record} !~ /\w+/) {
    $parser->next();
    $self->line_number($self->line_number + 1);
    $skipped++;
  }

  if($skipped) {
    my $ln = $self->line_number();

    $self->warning_msg(
      sprintf(
        "Skipped %i empty line%s from line %i",
        $skipped,
        $skipped > 1 ? 's' : '',
        $ln - $skipped + 1,
      )
    );
  }

  return $skipped;
}


=head2 line_number

  Arg 1      : (optional) int $line_number
  Example    : $ln = $parser->line_number();
  Description: Getter/setter for the current line number.
  Returntype : int
  Exceptions : none
  Caller     : next(), skip_empty_lines()
  Status     : Stable

=cut

sub line_number {
  my $self = shift;
  $self->{line_number} = shift if @_;
  return $self->{line_number};
}


=head2 valid_chromosomes

  Arg 1      : (optional) arrayref $valid_chromosomes
  Example    : $valids = $parser->valid_chromosomes();
  Description: Getter/setter for the list of valid chromosomes as found
               in the configured AnnotationSources.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : validate_vf()
  Status     : Stable

=cut

sub valid_chromosomes {
  my $self = shift;
  $self->{valid_chromosomes} = shift if @_;
  return $self->{valid_chromosomes};
}


=head2 delimiter

  Arg 1      : (optional) string $delimiter
  Example    : $delim = $parser->delimiter();
  Description: Getter/setter for the delimiter used to separate fields
               in the input. Most formats default to tab ("\t").
  Returntype : string
  Exceptions : none
  Caller     : detect_format(), parser()
  Status     : Stable

=cut

sub delimiter {
  my $self = shift;
  $self->{delimiter} = shift if @_;
  return $self->{delimiter};
}


=head2 detect_format

  Example    : $format = $parser->detect_format();
  Description: Attempts to detect the format of the input by analysing
               the first line. Not 100% successful, and will not work
               with STDIN as the filehandle cannot be "reset"
  Returntype : string
  Exceptions : throws if input is STDIN or cannot open file
  Caller     : new()
  Status     : Stable

=cut

sub detect_format {
  my $self = shift;

  my $file = $self->file;

  my $fh;
  my $file_was_fh = 0;

  if(openhandle($file)) {
    die("Cannot detect format from STDIN - specify format with --format [format]") if $file =~ /STDIN$/;
    $fh = $file;
    $file_was_fh = 1;
  }
  else {
    $fh = FileHandle->new();
    $fh->open($file) or throw("ERROR: Could not open $file to detect format\n");
  }

  my $format;
  my $delimiter = "\t";

  while(<$fh>) {
    next if /^\#/;
    chomp;

    # perl 5.8.8 doesn't recognise \v meaning any vertical whitespace...
    if($] < 5.01) {
      s/\r|(\x0D\x0A).+//;
    }
    else {
      s/\r|(?>\v|\x0D\x0A).+//;
    }

    # try to detect delimiter
    if(!/$delimiter/) {

      # try spaces
      if(/  /) {
        $delimiter = " +";
      }
      # try single space
      elsif(/ /) {
        $delimiter = " ";
      }
    }

    $self->delimiter($delimiter);
    $self->param('delimiter', $delimiter);

    my @data = split $delimiter, $_;
    next unless @data;

    # region chr21:10-10:1/A
    if (
      scalar @data == 1 &&
      $data[0] =~ /^[^\:]+\:\d+\-\d+(\:[\-\+]?1)?[\/\:](ins|dup|del|[ACGTN-]+)$/i
    ) {
      $format = 'region';
    }

    # SPDI: NC_000016.10:68684738:G:A
    elsif (
      scalar @data == 1 &&
      $data[0] =~ /^(.*?\:){2}([^\:]+|)$/i
    ) {
      $format = 'spdi';
    }

    # HGVS: ENST00000285667.3:c.1047_1048insC
    elsif (
      scalar @data == 1 &&
      $data[0] =~ /^([^\:]+)\:.*?([cgmrp]?)\.?([\*\-0-9]+.*)$/i
    ) {
      $format = 'hgvs';
    }

    # variant identifier: rs123456
    elsif (
      scalar @data == 1
    ) {
      $format = 'id';
    }

    # VCF: 20  14370  rs6054257  G  A  29  0  NS=58;DP=258;AF=0.786;DB;H2  GT:GQ:DP:HQ
    elsif (
      $data[0] =~ /(chr)?\w+/ &&
      $data[1] =~ /^\d+$/ &&
      $data[3] && $data[3] =~ /^[ACGTN\-\.]+$/i &&
      $data[4]
    ) {

      # do some more thorough checking on the ALTs
      my $ok = 1;

      foreach my $alt(split(',', $data[4])) {
        $ok = 0 unless $alt =~ /^[\.ACGTN\-\*]+$|^(\<[\w\:\*]+\>)$/i;
      }

      $format = 'vcf' if $ok;
    }

    # ensembl: 20  14370  14370  A/G  +
    elsif (
      $data[0] =~ /\w+/ &&
      $data[1] =~ /^\d+$/ &&
      $data[2] && $data[2] =~ /^\d+$/ &&
      $data[3] && $data[3] =~ /(ins|dup|del)|([ACGTN-]+\/[ACGTN-]+)/i
    ) {
      $format = 'ensembl';
    }

    # reset file handle if it was a handle
    eval {
      seek $fh, 0, 0 if $file_was_fh;
    };
    if($@) {
      close $fh;
      $self->file($self->{file_bak});
    }

    last;
  }

  return $format;
}


=head2 validate_vf

  Arg 1      : Bio::EnsEMBL::Variation::BaseVariationFeature
  Example    : $is_valid = $parser->validate_vf($vf);
  Description: Performs various (configurable) checks on a VariationFeature
               as produced by the parser:
                - creates a variation_name from the location+alleles if none
                - checks if chr is in user-specified list if provided
                - checks start and end look like numbers
                - checks start/end are valid (start <= end + 1) and corresponds to ref allele
                - checks chr is in valid list
                - checks ref allele vs genome if requested
  Returntype : bool
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub validate_vf {
  my $self = shift;
  my $vf = shift;

  # create name
  $vf->{variation_name} ||= sprintf(
    '%s_%i_%s',
    $vf->{original_chr} || $vf->{chr},
    $vf->{start},
    $vf->{allele_string} || $vf->{class_SO_term}
  );

  # user specified chr skip list
  if($self->{chr}) {
    return 0 unless grep {$vf->{chr} eq $_} @{$self->{chr}};
  }

  # sanity checks
  unless(looks_like_number($vf->{start}) && looks_like_number($vf->{end})) {
    $self->warning_msg("WARNING: Start ".$vf->{start}." or end ".$vf->{end}." coordinate invalid on line ".$self->line_number);
    return 0;
  }

  # check start <= end + 1
  if($vf->{start} > $vf->{end} + 1) {
    $self->warning_msg(
      "WARNING: start > end+1 : (START=".$vf->{start}.
      ", END=".$vf->{end}.
      ") on line ".$self->line_number."\n"
    );
    return 0;
  }

  # check we have this chr in any of the annotation sources
  # otherwise try to map to toplevel if available
  unless($self->_have_chr($vf)) {

    # slice adaptor required
    if(my $sa = $self->get_adaptor('core', 'Slice')) {
      $vf->{slice} ||= $self->get_slice($vf->{chr});

      if($vf->{slice}) {
        my $transformed = $vf->transform('toplevel');

        # copy to VF
        if($transformed) {
          $vf->{$_} = $transformed->{$_} for keys %$transformed;
          $vf->{original_chr} = $vf->{chr};
          $vf->{chr} = $vf->{slice}->seq_region_name;
        }

        # could not transform
        else {
          $self->warning_msg(
            "WARNING: Chromosome ".$vf->{chr}." not found in annotation sources or synonyms and could not transform to toplevel on line ".$self->line_number
          );
          return 0;
        }
      }

      # no slice
      else {
        $self->warning_msg(
          "WARNING: Could not fetch slice for chromosome ".$vf->{chr}." on line ".$self->line_number
        );
        return 0;
      }
    }

    # offline, can't transform
    else {
      $self->warning_msg(
        "WARNING: Chromosome ".$vf->{chr}." not found in annotation sources or synonyms on line ".$self->line_number
      );
      return 0;
    }
  }

  # structural variation?
  return $self->validate_svf($vf) if ref($vf) eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature';

  # uppercase allele string
  $vf->{allele_string} =~ tr/[a-z]/[A-Z]/;

  unless($vf->{allele_string} =~ /([ACGT-]+\/*)+/) {
    $self->warning_msg("WARNING: Invalid allele string ".$vf->{allele_string}." on line ".$self->line_number." or possible parsing error\n");
    return 0;
  }

  # insertion should have start = end + 1
  if($vf->{allele_string} =~ /^\-\// && $vf->{start} != $vf->{end} + 1) {
    my $variant_name = (defined $vf->name) ? "for variant (".$vf->name.") " : "";
    $self->warning_msg(
      "WARNING: Alleles look like an insertion (".
      $vf->{allele_string}.
      ") but coordinates are not start = end + 1 (START=".
      $vf->{start}.", END=".$vf->{end}.
      ") ".$variant_name."on line ".$self->line_number."\n"
    );
    return 0;
  }

  # check length of reference matches seq length spanned
  my @alleles = split '\/', $vf->{allele_string};

  # flag as unbalanced
  foreach my $allele(@alleles) {
    $allele =~ s/\-//g;
  }

  my $ref_allele = shift @alleles;

  if($ref_allele =~ /^[ACGT]*$/ && ($vf->{end} - $vf->{start}) + 1 != length($ref_allele)) {
    $self->warning_msg(
       "WARNING: Length of reference allele (".$ref_allele.
       " length ".length($ref_allele).") does not match co-ordinates ".$vf->{start}."-".$vf->{end}.
       " on line ".$self->line_number
    );
    return 0;
  }

  # check reference allele if requested
  if($self->{check_ref}) {
    my $ok = 0;
    my $slice_ref_allele;

    # insertion, therefore no ref allele to check
    if($ref_allele eq '') {
      $ok = 1;
    }
    else {
      my $slice_ref_allele = $self->_get_ref_allele($vf);

      if(!defined($slice_ref_allele)) {
        $self->warning_msg("WARNING: Could not fetch sub-slice from ".$vf->{chr}.":".$vf->{start}."\-".$vf->{end}."\(".$vf->{strand}."\) on line ".$self->line_number);
      }
      else {
        $ok = (uc($slice_ref_allele) eq uc($ref_allele) ? 1 : 0);
      }
    }

    if(!$ok) {
      $vf->{check_ref_failed} = 1;
      $self->warning_msg(
        "WARNING: Specified reference allele $ref_allele ".
        "does not match Ensembl reference allele".
        ($slice_ref_allele ? " $slice_ref_allele" : "").
        " on line ".$self->line_number
      );
      return 0;
    }
  }

  elsif($self->{lookup_ref}) {
    my $slice_ref_allele = $ref_allele eq '' ? '' : $self->_get_ref_allele($vf);

    if(!defined($slice_ref_allele)) {
      $self->warning_msg("WARNING: Could not reference allele for ".$vf->{chr}.":".$vf->{start}."\-".$vf->{end}."\(".$vf->{strand}."\) on line ".$self->line_number);
      return 0;
    }

    $vf->{allele_string} = join('/', map {$_ || '-'} ($slice_ref_allele, @alleles));
  }

  return 1;
}


=head2 _get_ref_allele

  Example    : $ref_allele = $parser->_get_ref_allele($vf);
  Description: Looks up the reference sequence covered by the given
               VariationFeature. Uses a slice that will fetch from
               FASTA (if configured) or database (assuming not in
               offline mode). Returns undef if it failed to fetch an
               appropriate slice.
  Returntype : string
  Exceptions : none
  Caller     : validate_vf()
  Status     : Stable

=cut

sub _get_ref_allele {
  my ($self, $vf) = @_;

  $vf->{slice} ||= $self->get_slice($vf->{chr});

  my $slice_ref = $vf->{slice}->sub_Slice($vf->{start}, $vf->{end}, $vf->{strand});

  return $slice_ref ? $slice_ref->seq : undef;
}


=head2 _have_chr

  Example    : $have_chr = $parser->_have_chr();
  Description: Checks if the chromosome name assigned in the VariationFeature
               is valid given the list provided by valid_chromosomes(); also
               attempts to transform by adding/removing "chr".
  Returntype : bool
  Exceptions : none
  Caller     : validate_vf()
  Status     : Stable

=cut

sub _have_chr {
  my ($self, $vf) = @_;

  my $have_chr = 0;
  my $valid = $self->valid_chromosomes;
  my $vf_chr = $vf->{chr};

  if($valid->{$vf_chr}) {
    $have_chr = 1;
  }
  else {
    my $synonyms = $self->chromosome_synonyms;
    $valid->{$_} = 1 for map {keys %{$synonyms->{$_} || {}}} keys %$valid;

    if($valid->{$vf->{chr}}) {
      $have_chr = 1;
    }
    else {
      my $tmp_chr = $vf_chr;
      $tmp_chr =~ s/^chr//i;

      if($valid->{'chr'.$vf->{chr}} || $valid->{$tmp_chr}) {
        $have_chr = 1;
      }
      elsif($vf->{chr} eq 'M' && $valid->{MT}) {
        $vf->{chr} = 'MT';
        $have_chr = 1;
      }
    }
  }

  return $have_chr;
}


=head2 validate_svf

  Arg 1      : Bio::EnsEMBL::Variation::StructuralVariationFeature
  Example    : $is_valid = $parser->validate_svf($svf);
  Description: Stub, not currently implemented
  Returntype : bool
  Exceptions : none
  Caller     : validate_vf()
  Status     : Stable

=cut

sub validate_svf {
  return 1;
}


=head2 post_process_vfs

  Arg 1      : arrayref of Bio::EnsEMBL::Variation::BaseVariationFeature
  Example    : $vfs = $parser->post_process_vfs($svf);
  Description: Applies some optional post-processing to VariationFeatures:
                - mapping to LRG (--lrg)
                - minimising alleles (--minimal)
  Returntype : arrayref of Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub post_process_vfs {
  my $self = shift;
  my $vfs = shift;

  # map to LRGs
  $vfs = $self->map_to_lrg($vfs) if $self->{lrg};

  # minimise alleles?
  $vfs = $self->minimise_alleles($vfs) if $self->{minimal};

  # copy start, end coords to seq_region_start, seq_region_end
  # otherwise for circular chromosomes the core API will try to do a DB lookup and die
  foreach my $vf(@$vfs) {
    $vf->seq_region_start($vf->{start});
    $vf->seq_region_end($vf->{end});
  }

  return $vfs;
}


=head2 map_to_lrg

  Arg 1      : arrayref of Bio::EnsEMBL::Variation::BaseVariationFeature
  Example    : $vfs = $parser->map_to_lrg($svf);
  Description: Maps variants to LRGs, appending any successfully mapped
               VariationFeatures to the input arrayref.
  Returntype : arrayref of Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : none
  Caller     : post_process_vfs()
  Status     : Stable

=cut

sub map_to_lrg {
  my $self = shift;
  my $vfs = shift;

  return $vfs if $self->{dont_map_to_lrg};

  my @return;

  foreach my $vf(@$vfs) {

    # add the unmapped VF to the array
    push @return, $vf;

    # make sure the VF has an attached slice
    $vf->{slice} ||= $self->get_slice($vf->{chr});
    unless(defined($vf->{slice})) {
      push @return, $vfs;
      next;
    }

    if($self->{can_map_to_lrg} || $vf->{slice}->coord_system->adaptor->fetch_by_name('lrg')) {
      $self->{can_map_to_lrg} = 1;
    }
    else {
      $self->{dont_map_to_lrg} = 1;
      return $vfs;
    }

    # transform LRG <-> chromosome
    my $new_vf;

    eval { $new_vf = $vf->transform($vf->{slice}->coord_system->name eq 'lrg' ? 'chromosome' : 'lrg') };

    # add it to the array if transformation worked
    if(defined($new_vf)) {

      # update new VF's chr entry
      $new_vf->{chr} = $new_vf->seq_region_name;
      push @return, $new_vf;
    }
  }

  return \@return;
}


=head2 minimise_alleles

  Arg 1      : arrayref of Bio::EnsEMBL::Variation::BaseVariationFeature
  Example    : $vfs = $parser->minimise_alleles($svf);
  Description: Modifies VariationFeatures by reducing REF/ALT to their
               minimal shared sequence. Same principle as
               InputBuffer::split_variants, but runs on VFs with a single ALT.

               The original allele_string, start and end keys are stored
               on the modified VF as e.g. original_allele_string.
  Returntype : arrayref of Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : none
  Caller     : post_process_vfs()
  Status     : Stable

=cut

sub minimise_alleles {
  my $self = shift;
  my $vfs = shift;

  my @return;

  foreach my $vf(@$vfs) {

    # skip VFs with more than one alt
    # they get taken care of later by split_variants/rejoin_variants
    if(!$vf->{allele_string} || $vf->{allele_string} =~ /.+\/.+\/.+/ || $vf->{allele_string} !~ /.+\/.+/) {
      push @return, $vf;
    }

    else {
      my @alleles = split('/', $vf->{allele_string});
      my $ref = shift @alleles;

      foreach my $alt(@alleles) {

        my $start = $vf->{start};
        my $end   = $vf->{end};

        ($ref, $alt, $start, $end) = @{trim_sequences($ref, $alt, $start, $end, 1)};

        # create a copy
        my $new_vf;
        %$new_vf = %{$vf};
        bless $new_vf, ref($vf);

        # give it a new allele string and coords
        $new_vf->allele_string($ref.'/'.$alt);
        $new_vf->{start}                  = $start;
        $new_vf->{end}                    = $end;
        $new_vf->{seq_region_start}       = $start;
        $new_vf->{seq_region_end}         = $end;
        $new_vf->{original_allele_string} = $vf->{allele_string};
        $new_vf->{original_start}         = $vf->{start};
        $new_vf->{original_end}           = $vf->{end};
        $new_vf->{minimised}              = 1;

        push @return, $new_vf;
      }
    }
  }

  return \@return;
}

1;
