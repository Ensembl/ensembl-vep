=head1 LICENSE

Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::VCF - VCF input parser

=head1 SYNOPSIS

my $parser = Bio::EnsEMBL::VEP::Parser::VCF->new({
  config => $config,
  file   => 'variant.vcf',
});

my $vf = $parser->next();

=head1 DESCRIPTION

VCF format parser.

4.2 spec: https://samtools.github.io/hts-specs/VCFv4.2.pdf

IMPORTANT NOTE: unbalanced substitutions are encoded in VCF
with the preceding base prepended to the REF and ALT alleles.
In order to be processed by VEP and the Ensembl API, these are
converted to Ensembl style by removing the prepended base and
incrementing the start position, but only if the same first base
is shared by *all* REF and ALT alleles.

This behaviour can be modified further by using the --minimal
flag, which reduces each REF+ALT pair to their minimal shared
sequence.

=head1 METHODS

=cut


use strict;
use warnings;
no warnings 'recursion';

package Bio::EnsEMBL::VEP::Parser::VCF;

use base qw(Bio::EnsEMBL::VEP::Parser);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::Parser::VCF4;

use List::Util qw(max);
use List::MoreUtils qw(uniq);

=head2 new

  Arg 1      : hashref $args
               {
                 config    => Bio::EnsEMBL::VEP::Config,
                 file      => string or filehandle,
               }
  Example    : $parser = Bio::EnsEMBL::VEP::Parser::VCF->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::Parser::VCF object.
  Returntype : Bio::EnsEMBL::VEP::Parser::VCF
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(allow_non_variant gp individual individual_zyg process_ref_homs phased)]);

  return $self;
}


=head2 parser

  Example    : $io_parser = $parser->parser();
  Description: Get ensembl-io parser object used to read data from input.
  Returntype : Bio::EnsEMBL::IO::Parser::VCF4
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;

  if(!exists($self->{parser})) {
    $self->{parser} = Bio::EnsEMBL::IO::Parser::VCF4->open($self->file);
    $self->{parser}->{delimiter} = $self->delimiter;
  }

  return $self->{parser};
}


=head2 headers

  Example    : $headers = $parser->headers();
  Description: Gets headers from the VCF input. Required for
               writing output as VCF.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub headers {
  my $self = shift;

  if(!exists($self->{headers})) {
    my $parser = $self->parser;

    unless($self->{_have_read_next}) {
      $parser->next;
      $self->{_have_read_next} = 1;
    }

    $self->{headers} = $parser->{_raw_metadata} || [];
  }

  return $self->{headers};
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

    # getting the header requires we trigger next once
    # so we don't want to trigger it again (once)
    if($self->{_have_read_next}) {
      delete $self->{_have_read_next};
    }
    else {
      $self->parser->next;
    }

    push @$cache, @{$self->create_VariationFeatures()};
  }

  my $vf = shift @$cache;
  return $vf unless $vf;

  unless($self->validate_vf($vf) || $self->{dont_skip}) {
    return $self->next();
  }

  return $vf;
}


=head2 create_VariationFeatures

  Example    : $vfs = $parser->create_VariationFeatures();
  Description: Create one or more VariationFeature objects from the current line
               of input; multiple may be returned if multiple individuals
               are requested using --individual.
  Returntype : arrayref of Bio::EnsEMBL::BaseVariationFeature
  Exceptions : warns if GP flag required and not found
  Caller     : next()
  Status     : Stable

=cut

sub create_VariationFeatures {
  my $self = shift;

  my $parser = $self->parser();
  my $record = $parser->{record};

  return [] unless $record && @$record;

  $self->line_number($self->line_number + 1);

  # get the data we need to decide if this is an SV
  my ($ref, $alts, $info) = (
    $parser->get_reference,
    $parser->get_alternatives,
    $parser->get_info,
  );

  if(
    join("", $ref, @$alts) !~ /^[ACGT]+$/ &&
    (
      $info->{SVTYPE} || # deprecated in VCF 4.4
      join(",", @$alts) =~ /[<\[\]][^\*]+[>\]\[]|^\.\w+|\w+\.$/
    )
  ) {
    return $self->create_StructuralVariationFeatures();
  }

  # get the rest of the relevant data
  my ($chr, $start, $end, $ids) = (
    $parser->get_seqname,
    $parser->get_raw_start,
    $parser->get_raw_start,
    $parser->get_IDs,
  );
  $end += length($ref) - 1;

  # non-variant
  my $non_variant = 0;

  if($alts->[0] eq '.') {
    if($self->{allow_non_variant}) {
      $non_variant = 1;
    }
    else {
      $parser->next();
      return $self->create_VariationFeatures;
    }
  }

  # some VCF files have a GRCh37 pos defined in GP flag in INFO column
  # if user has requested, we can use that as the position instead
  if($self->{gp}) {
    $chr = undef;
    $start = undef;

    if(my $gp = $info->{GP}) {
      ($chr, $start) = split ':', $gp;
      $end = $start;
    }

    unless(defined($chr) and defined($start)) {
      $self->skipped_variant_msg("No GP flag found in INFO column");
      $parser->next();
      return $self->create_VariationFeatures;
    }
  }

  # record original alleles
  # if they get changed, we need to map from old to new in create_individual_VariationFeatures
  my @original_alleles = ($ref, @$alts);

  # adjust end coord
  # $end += (length($ref) - 1);

  # find out if any of the alt alleles make this an insertion or a deletion
  my $is_indel = 0;
  foreach my $alt_allele(@$alts) {
    $is_indel = 1 if $alt_allele =~ /^[DI]/ or length($alt_allele) != length($ref);
  }

  # multiple alt alleles?
  if(scalar @$alts > 1) {
    if($is_indel) {

      # find out if all the alts start with the same base
      # ignore "*"-types
      my %first_bases = map {substr($_, 0, 1) => 1} grep {!/\*/} ($ref, @$alts);

      if(scalar keys %first_bases == 1) {
        $ref = substr($ref, 1) || '-';
        $start++;

        my @new_alts;

        foreach my $alt_allele(@$alts) {
          $alt_allele = substr($alt_allele, 1) unless $alt_allele =~ /\*/;
          $alt_allele = '-' if $alt_allele eq '';
          push @new_alts, $alt_allele;
        }

        $alts = \@new_alts;
      }
    }
  }

  elsif($is_indel) {

    # insertion or deletion (VCF 4+)
    if(substr($ref, 0, 1) eq substr($alts->[0], 0, 1)) {

      # chop off first base
      $ref = substr($ref, 1) || '-';
      $alts->[0] = substr($alts->[0], 1) || '-';

      $start++;
    }
  }

  # create VF object
  my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => $start,
    end            => $end,
    allele_string  => $non_variant ? $ref : $ref.'/'.join('/', @$alts),
    strand         => 1,
    map_weight     => 1,
    adaptor        => $self->get_adaptor('variation', 'VariationFeature'),
    variation_name => @$ids ? $ids->[0] : undef,
    chr            => $chr,
    _line          => $record,
  });

  # flag as non-variant
  $vf->{non_variant} = 1 if $non_variant;

  # individual data?
  if($self->{individual} || $self->{individual_zyg}) {
    my @changed_alleles = ($ref, @$alts);
    my %allele_map = map {$original_alleles[$_] => $changed_alleles[$_]} 0..$#original_alleles;
    
    my $method = $self->{individual} ? 'create_individual_VariationFeatures' : 'create_individuals_zyg_VariationFeature';

    my @return =
      map {@{$self->$method($_, \%allele_map)}}
      @{$self->post_process_vfs([$vf])};

    # if all selected individuals had REF or missing genotypes @return will be empty
    # re-run this method in this case to avoid the parser shorting out and finishing
    if(@return) {
      return \@return;
    }
    else {
      $parser->next();
      return $self->create_VariationFeatures;
    }
  }

  # normal return
  return $self->post_process_vfs([$vf]);
}


sub _expand_tandem_repeat_allele_string {
  my ($self, $type, $info, $ref) = @_;

  if ($info->{RUS} && ($info->{RUC} || $info->{RB}) ) {
    # warn that CIRUC and CIRB INFO fields are ignored
    $self->warning_msg("CIRUC and CIRB INFO fields are ignored when calculating alternative alleles in tandem repeats")
      if $info->{CIRUC} or $info->{CIRB};

    # RN  : Total number of repeat sequences in this allele
    # RUS : Repeat unit sequence of the corresponding repeat sequence
    # RUC : Repeat unit count of corresponding repeat sequence
    # RB  : Total number of bases in the corresponding repeat sequence
    my @RN = $info->{RN} ?
      split(/,/, $info->{RN}) :
      # if RN is missing, assume 1 repeat for each CNV:TR allele
      (1) x ( () = $type =~ /CNV:TR/g );
    my @RUS = split /,/,  $info->{RUS};
    my @RUC = split /,/, ($info->{RUC} || ''); # undefined if no RUC
    my @RB  = split /,/, ($info->{RB}  || ''); # undefined if no RB
    my $is_RUC_defined = scalar @RUC;

    my $max_sv_size = $self->param('max_sv_size') || 5000;
    my $is_oversized = 0;

    # build sequence of tandem repeat based on INFO fields
    my @repeats;
    for my $records (@RN) {
      my $repeat;
      next if $records eq '.'; # ignore missing values
      for ( 1 .. $records ) {
        my $seq  = shift(@RUS);
           $seq  = 'N' if $seq eq '.'; # missing sequence
        my $num  = $is_RUC_defined ? shift(@RUC) : shift(@RB) / length($seq);

        # avoid calculating really large repeats
        $is_oversized = length($seq) * $num > $max_sv_size;
        last if $is_oversized;
        $repeat .= $seq x $num;
      }

      last if $is_oversized;
      # prepend reference allele
      push @repeats, $repeat;
    }

    # avoid storing alternative allele for tandem repeats with large repeats
    return $ref.'/'.join('/', @repeats) unless $is_oversized;
  }
  return undef;
}


=head2 create_StructuralVariationFeatures

  Example    : $vfs = $parser->create_StructuralVariationFeatures();
  Description: Create StructuralVariationFeature objects from the current line
               of input.
  Returntype : arrayref of Bio::EnsEMBL::StructuralVariationFeature
  Exceptions : warns if SV type indicated but end coord can't be determined
  Caller     : create_VariationFeatures()
  Status     : Stable

=cut

sub create_StructuralVariationFeatures {
  my $self = shift;

  my $parser = $self->parser();
  my $record = $parser->{record};
  my $skip_line;

  # get relevant data
  my ($chr, $start, $end, $ref, $alts, $info, $ids) = (
    $parser->get_seqname,
    $parser->get_start,
    $parser->get_end,
    $parser->get_reference,
    $parser->get_alternatives,
    $parser->get_info,
    $parser->get_IDs,
  );

  ## get structural variant type from ALT or (deprecated) SVTYPE tag
  my $alt = join("/", @$alts);
  my $type = $alt;

  # replace with SVTYPE tag if ALT does not follow VCF 4.4 specs
  if ($info->{SVTYPE} && $alt !~ /^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)/) {
    $type = $info->{SVTYPE};
  }

  my $so_term = $self->get_SO_term($type);
  unless ($so_term) {
    $skip_line = 1;
    $so_term   = $type;
  }

  ## get breakends from INFO field (from Illumina Manta, for instance)
  if ($so_term =~ /breakpoint/) {
    ## Illumina Manta (SV caller) may use INFO/END to identify the position of
    ## the breakend mate (this is not supported by VCF 4.4 specifications)
    my $incorrect_end = $info->{END};
    delete $parser->get_info->{END};
    $end = $parser->get_end if $incorrect_end;

    if (defined $info->{CHR2}) {
      my $breakend_chr = $self->get_source_chr_name($info->{CHR2});
      my $breakend_pos = $info->{END2} || $incorrect_end;
      if (defined $breakend_chr and defined $breakend_pos) {
        $alt = $alt =~ /^<?BND>?$/i ? "N" : "$alt/N";
        $alt = sprintf('%s[%s:%s[', $alt, $breakend_chr, $breakend_pos);
      }
    }
    $alt = $ref . "/$alt" unless $alt =~ /^\.|\.$/;
  }
  ## parse tandem repeats based on VCF INFO fields
  elsif ($so_term =~ /tandem/ ) {
    # validate SVLEN
    if (defined $info->{SVLEN}) {
      my @svlen = uniq(split(/,/, $info->{SVLEN}));
      # warn if there are different references per alternative allele
      $self->warning_msg(
        "found tandem repeats with different references per alternative allele: " .
        "SVLEN=" . $info->{SVLEN} . "; only using reference with largest size"
      ) if scalar @svlen > 1;
      $end = $start + max(@svlen) - 1;
    }

    # get reference allele and allele_string from tandem repeats
    my $allele_string;
    my $slice_ref = $self->get_slice($chr);
    if ($slice_ref) {
      $slice_ref = $slice_ref->sub_Slice($start, $end, 1);
      $ref = $slice_ref->seq if defined $slice_ref and defined $slice_ref->seq;
      $allele_string = $self->_expand_tandem_repeat_allele_string($type, $info, $ref) if defined $ref;
    } elsif (!defined($self->param('fasta')) && $self->param('offline')) {
      $self->warning_msg(
        "could not fetch sequence for tandem repeats; " .
        "consequence calculation for tandem repeats is less precise when using --offline without --fasta\n");
    } else {
      $self->warning_msg("could not fetch sequence for tandem repeat");
    }

    if (defined $allele_string) {
      # convert tandem repeats to VariationFeature in order to properly
      # calculate consequences based on alternative allele sequence
      my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
        start          => $start,
        end            => $end,
        allele_string  => $allele_string,
        strand         => 1,
        map_weight     => 1,
        adaptor        => $self->get_adaptor('variation', 'VariationFeature'),
        variation_name => @$ids ? $ids->[0] : undef,
        chr            => $chr,
        _line          => $record,
      });
      return $self->post_process_vfs([$vf]);
    }
  }

  # check for imprecise breakpoints
  my ($min_start, $max_start, $min_end, $max_end) = (
    $parser->get_outer_start,
    $parser->get_inner_start,
    $parser->get_inner_end,
    $parser->get_outer_end,
  );

  if($start >= $end && $so_term =~ /del/i) {
    $self->skipped_variant_msg("deletion looks incomplete");
    $skip_line = 1;
  }

  my $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
    start          => $start,
    inner_start    => $max_start,
    outer_start    => $min_start,
    end            => $end,
    inner_end      => $min_end,
    outer_end      => $max_end,
    strand         => 1,
    adaptor        => $self->get_adaptor('variation', 'StructuralVariationFeature'),
    variation_name => @$ids ? $ids->[0] : undef,
    chr            => $self->get_source_chr_name($chr),
    class_SO_term  => $so_term,
    allele_string  => $alt,
    _line          => $record
  });
  $svf->{vep_skip} = $skip_line if defined $skip_line;
  return $self->post_process_vfs([$svf]);
}


=head2 create_individual_VariationFeatures
 
  Arg 1      : Bio::EnsEMBL::VariationFeature $vf
  Arg 2      : hashref $allele_map
  Example    : $vfs = $parser->create_individual_VariationFeatures($vf, $map);
  Description: Create one VariationFeature object per configured
               individual/sample. Arg 2 $allele_map is a hashref mapping the
               allele index to the actual ALT string it represents.
  Returntype : arrayref of Bio::EnsEMBL::VariationFeature
  Exceptions : none
  Caller     : create_VariationFeatures()
  Status     : Stable

=cut

sub create_individual_VariationFeatures {
  my $self = shift;
  my $vf = shift;
  my $allele_map = shift;

  my $parser = $self->parser();
  my $record = $parser->{record};

  my @alleles = split '\/', $vf->{allele_string};
  my $ref = $alleles[0];

  my @return;

  # get genotypes from parser
  my $include = lc($self->{individual}->[0]) eq 'all' ? $parser->get_samples : $self->{individual};

  # Compare sample names
  if(lc($self->{individual}->[0]) ne 'all') {
    my $found = _find_in_array($parser->get_samples, $include);

    if(!$found) {
      die("ERROR: Sample IDs given (", join(",", @{$include}), ") do not match samples from VCF (", join(",", @{$parser->get_samples}), ")\n");
    }
  }

  my $ind_gts = $parser->get_samples_genotypes($include, 1 - ($self->{allow_non_variant} || 0));

  foreach my $ind(@$include) {

    # get alleles present in this individual
    my $gt = $ind_gts->{$ind};
    next if (!$gt);
    my @bits = map { $allele_map->{$_} } split /\||\/|\\/, $gt;
    my $phased = ($gt =~ /\|/ ? 1 : 0);

    # shallow copy VF
    my $vf_copy = { %$vf };
    bless $vf_copy, ref($vf);

    # get non-refs, remembering to exclude "*"-types
    my %non_ref = map {$_ => 1} grep {$_ ne $ref && $_ !~ /\*/} @bits;

    # construct allele_string
    if(scalar keys %non_ref) {
      $vf_copy->{allele_string} = $ref."/".(join "/", keys %non_ref);
    }
    else {
      $vf_copy->{allele_string} = $ref;
      $vf_copy->{hom_ref} = 1;

      if($self->{process_ref_homs}) {
        $vf_copy->{allele_string} .= "/".$ref;
      }
      else {
        $vf_copy->{non_variant} = 1;
      }
    }

    # store phasing info
    $vf_copy->{phased} = $self->{phased} ? 1 : $phased;

    # store GT
    $vf_copy->{genotype} = \@bits;

    # store individual name
    $vf_copy->{individual} = $ind;

    push @return, $vf_copy;
  }

  return \@return;
}

=head2 create_individuals_zyg_VariationFeature
 
  Arg 1      : Bio::EnsEMBL::VariationFeature $vf
  Arg 2      : hashref $allele_map
  Example    : $vfs = $parser->create_individuals_zyg_VariationFeature($vf, $map);
  Description: Create one VariationFeature object with
               individual/sample info. Arg 2 $allele_map is a hashref mapping the
               allele index to the actual ALT string it represents.
  Returntype : arrayref of Bio::EnsEMBL::VariationFeature
  Exceptions : none
  Caller     : create_VariationFeatures()
  Status     : Stable

=cut

sub create_individuals_zyg_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $allele_map = shift;

  my $parser = $self->parser();
  my $record = $parser->{record};

  my @alleles = split '\/', $vf->{allele_string};
  my $ref = $alleles[0];

  my @return;

  # get genotypes from parser
  my $include = lc($self->{individual_zyg}->[0]) eq 'all' ? $parser->get_samples : $self->{individual_zyg};

  # Compare sample names
  if(lc($self->{individual_zyg}->[0]) ne 'all') {
    my $found = _find_in_array($parser->get_samples, $include);

    if(!$found) {
      die("ERROR: Sample IDs given (", join(",", @{$include}), ") do not match samples from VCF (", join(",", @{$parser->get_samples}), ")\n");
    }
  }

  my $ind_gts = $parser->get_samples_genotypes($include, 1 - ($self->{allow_non_variant} || 0));

  my $n_individuals = scalar(@{$include});

  foreach my $ind(@$include) {
    # get alleles present in this individual
    my $gt = $ind_gts->{$ind};
    next if (!$gt);
    my @bits = map { $allele_map->{$_} } split /\||\/|\\/, $gt;
    my $phased = ($gt =~ /\|/ ? 1 : 0);

    # get non-refs, remembering to exclude "*"-types
    my %non_ref = map {$_ => 1} grep {$_ ne $ref && $_ !~ /\*/} @bits;

    # Genotype is reference
    if(!scalar keys %non_ref) {
      $vf->{hom_ref}->{$ind} = 1;
      $vf->{non_variant}->{$ind} = 1;

      if($n_individuals == 1) {
        $vf->{allele_string} = $ref."/".$ref ;
        $vf->{unique_ind} = 1;
      }
    }

    # store phasing info
    $vf->{phased}->{$ind} = $self->{phased} ? 1 : $phased;

    # store GT
    $vf->{genotype_ind}->{$ind} = \@bits;
  }

  push @return, $vf;

  return \@return;
}

sub _find_in_array {
  my $all_samples = shift;
  my $samples = shift;

  my %all = map { $_ => 1 } @{$all_samples};
  my %input_samples = map { $_ => 1 } @{$samples};

  foreach my $sample (keys %input_samples) {
    if(!$all{$sample}) {
      return 0;
    }
  }

  return 1;
}

1;
