=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Parser;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(trim_sequences get_compressed_filehandle);

use Bio::EnsEMBL::VEP::Parser::VCF;
use Bio::EnsEMBL::VEP::Parser::VEP_input;
use Bio::EnsEMBL::VEP::Parser::ID;
use Bio::EnsEMBL::VEP::Parser::HGVS;

use Scalar::Util qw(openhandle looks_like_number);
use FileHandle;

my %FORMAT_MAP = (
  'vcf'     => 'VCF',
  'ensembl' => 'VEP_input',
  'id'      => 'ID',
  'hgvs'    => 'HGVS',
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  $self->add_shortcuts([qw(dont_skip check_ref chr lrg minimal delimiter)]);

  my $hashref = $_[0];

  throw("ERROR: No file given\n") unless $hashref->{file};
  $self->file($hashref->{file});

  $self->delimiter($hashref->{delimiter}) if $hashref->{delimiter};

  $self->line_number(0);

  $self->valid_chromosomes({map {$_ => 1} @{$hashref->{valid_chromosomes} || []}});

  if(my $format = $hashref->{format}) {

    delete $hashref->{format};

    # detect format
    if(lc($format) eq 'guess' || lc($format) eq 'detect' || lc($format) eq 'auto') {
      $format = $self->detect_format();
    }

    throw("ERROR: Can't detect format\n") unless $format;

    $format = lc($format);
    throw("ERROR: Unknown or unsupported format $format\n") unless $FORMAT_MAP{$format};

    $self->param('format', $format);

    my $class = 'Bio::EnsEMBL::VEP::Parser::'.$FORMAT_MAP{$format};
    return $class->new({%$hashref, config => $self->config, delimiter => $self->delimiter});
  }

  return $self;
}

# generic next method
# sub-classes may override it
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

sub headers {
  return [];
}

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
      $file = get_compressed_filehandle($file, 1) if -B $file;
    }

    $self->{file} = $file;
  }
  
  return $self->{file};
}

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

sub line_number {
  my $self = shift;
  $self->{line_number} = shift if @_;
  return $self->{line_number};
}

sub valid_chromosomes {
  my $self = shift;
  $self->{valid_chromosomes} = shift if @_;
  return $self->{valid_chromosomes};
}

sub delimiter {
  my $self = shift;
  $self->{delimiter} = shift if @_;
  return $self->{delimiter};
}

sub detect_format {
  my $self = shift;

  my $file = $self->file;

  my $fh;
  my $file_was_fh = 0;

  if(openhandle($file)) {
    throw("Cannot detect format from STDIN - specify format with --format [format]") if $file =~ /STDIN$/;
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
    s/\r|(?>\v|\x0D\x0A).+//;

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

    # HGVS: ENST00000285667.3:c.1047_1048insC
    if (
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

    # pileup: chr1  60  T  A
    elsif (
      $data[0] =~ /(chr)?\w+/ &&
      $data[1] =~ /^\d+$/ &&
      $data[2] && $data[2] =~ /^[\*ACGTN-]+$/i &&
      $data[3] && $data[3] =~ /^[\*ACGTNRYSWKM\+\/-]+$/i
    ) {
      $format = 'pileup';
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

# takes VFs created from input, fixes and checks various things
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
    $self->warning_msg(
      "WARNING: Alleles look like an insertion (".
      $vf->{allele_string}.
      ") but coordinates are not start = end + 1 (START=".
      $vf->{start}.", END=".$vf->{end}.
      ") on line ".$self->line_number."\n"
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
      $vf->{slice} ||= $self->get_slice($vf->{chr});

      my $slice_ref = $vf->{slice}->sub_Slice($vf->{start}, $vf->{end}, $vf->{strand});

      if(!defined($slice_ref)) {
        $self->warning_msg("WARNING: Could not fetch sub-slice from ".$vf->{chr}.":".$vf->{start}."\-".$vf->{end}."\(".$vf->{strand}."\) on line ".$self->line_number);
      }

      else {
        $slice_ref_allele = $slice_ref->seq;
        $ok = (uc($slice_ref_allele) eq uc($ref_allele) ? 1 : 0);
      }
    }

    if(!$ok) {
      $self->warning_msg(
        "WARNING: Specified reference allele $ref_allele ".
        "does not match Ensembl reference allele".
        ($slice_ref_allele ? " $slice_ref_allele" : "").
        " on line ".$self->line_number
      );
      return 0;
    }
  }

  return 1;
}

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

# validate a structural variation
sub validate_svf {
  return 1;
}

# post process VFs - remapping etc
sub post_process_vfs {
  my $self = shift;
  my $vfs = shift;

  # map to LRGs
  $vfs = $self->map_to_lrg($vfs) if $self->{lrg};

  # minimise alleles?
  $vfs = $self->minimise_alleles($vfs) if $self->{minimal};

  return $vfs;
}

sub map_to_lrg {
  my $self = shift;
  my $vfs = shift;

  my @return;

  foreach my $vf(@$vfs) {

    # add the unmapped VF to the array
    push @return, $vf;

    # make sure the VF has an attached slice
    $vf->{slice} ||= $self->get_slice($vf->{chr});
    next unless defined($vf->{slice});

    # transform LRG <-> chromosome
    my $new_vf = $vf->transform($vf->{slice}->coord_system->name eq 'lrg' ? 'chromosome' : 'lrg');

    # add it to the array if transformation worked
    if(defined($new_vf)) {

      # update new VF's chr entry
      $new_vf->{chr} = $new_vf->seq_region_name;
      push @return, $new_vf;
    }
  }

  return \@return;
}

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

        ($ref, $alt, $start, $end) = @{trim_sequences($ref, $alt, $start, $end)};
        $ref ||= '-';
        $alt ||= '-';

        # create a copy
        my $new_vf;
        %$new_vf = %{$vf};
        bless $new_vf, ref($vf);

        # give it a new allele string and coords
        $new_vf->allele_string($ref.'/'.$alt);
        $new_vf->{start}                  = $start;
        $new_vf->{end}                    = $end;
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
