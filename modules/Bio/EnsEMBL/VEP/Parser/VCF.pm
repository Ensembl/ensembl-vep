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
      join(",", @$alts) =~ /[<\[\]][^\*]+[>\]\[]/
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


=head2 _parse_breakend_alt_alleles

  Arg 1      : string $alt (ALT field from VCF)
  Arg 2      : hashref $info (INFO field from VCF)
  Example    : $parsed_alleles = $self->_parse_breakend_alt_alleles($alt, $info);
  Description: Parse alternative alleles from breakend structural variants
  Returntype : Arrayref of hashrefs containing information of the parsed alleles
  Exceptions : none
  Caller     : create_StructuralVariationFeatures()
  Status     : Stable

=cut

sub _parse_breakend_alt_alleles {
  my $self = shift;
  my $alt = shift;
  my $info = shift;

  my $parsed_allele = [];

  # Support multiple breakends from ALT in the format C[2:321682[,]17:198982]G
  for my $alt_string (split ",", $alt) {
    my ($alt_allele) = ($alt_string =~ '[\[\]]?([A-Za-z]+)[\[\]]?');
    my ($alt_chr, $alt_pos) = ($alt_string =~ '([A-Za-z0-9]+) ?: ?([0-9]+)');

    next if !defined $alt_allele or !defined $alt_chr or !defined $alt_pos;
    $alt_chr = $self->get_source_chr_name($alt_chr);

    # Mate orientation:
    #   C[2:321682[  - normal
    #   G]17:198982] - inverted
    #   ]13:123456]T - normal
    #   [17:198983[A - inverted
    my $inverted = ($alt_string =~ '^\[|\]$') || 0;

    # Mate placement:
    #   C[2:321682[  - left
    #   [17:198983[A - right
    my $placement = ($alt_string =~ '^(\[|\])') ? 'left' : 'right';

    my $slice = Bio::EnsEMBL::Slice->new_fast({
      seq_region_name => $alt_chr,
      start => $alt_pos,
      end => $alt_pos,
    });

    push @$parsed_allele, {
      string    => $alt_string,
      allele    => $alt_allele,
      pos       => $alt_pos,
      start     => $alt_pos,
      end       => $alt_pos,
      chr       => $alt_chr,
      inverted  => $inverted,
      placement => $placement,
      slice     => $slice
    };
  }

  # Support breakend described in INFO field
  if (defined $info->{CHR2} && scalar @$parsed_allele == 0) {
    my $alt_chr = $self->get_source_chr_name($info->{CHR2});
    my $alt_pos = $info->{END2};

    if (defined $alt_chr && defined $alt_pos) {
      my $slice = Bio::EnsEMBL::Slice->new_fast({
        seq_region_name => $alt_chr,
        start => $alt_pos,
        end => $alt_pos,
      });

      push @$parsed_allele, {
        string    => $alt,
        allele    => $alt,
        pos       => $alt_pos,
        start     => $alt_pos,
        end       => $alt_pos,
        chr       => $alt_chr,
        slice     => $slice
      };
    }
  }

  return $parsed_allele;
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
  my ($chr, $start, $end, $alts, $info, $ids) = (
    $parser->get_seqname,
    $parser->get_start,
    $parser->get_end,
    $parser->get_alternatives,
    $parser->get_info,
    $parser->get_IDs,
  );
  
  ## get structural variant type from SVTYPE tag (deprecated in VCF 4.4) or ALT
  my $alt = join(",", @$alts);
  my $type = $alt ne '.' ? $alt : $info->{SVTYPE};
  my $so_term = $self->get_SO_term($type);

  unless ($so_term) {
    $skip_line = 1;
    $so_term   = $type;
  }

  ## Illumina Manta (SV caller) may use INFO/END to identify the position of the
  ## breakend mate (unsupported based on VCF 4.4 specifications)
  if ($so_term =~ /break/ && $parser->get_info->{END}) {
    # INFO/END2 is used to identify the position of the breakend mate
    $parser->get_info->{END2} ||= $parser->get_info->{END};
    delete $parser->get_info->{END};
    $end = $parser->get_end;
  }

  # check for imprecise breakpoints
  my ($min_start, $max_start, $min_end, $max_end) = (
    $parser->get_outer_start,
    $parser->get_inner_start,
    $parser->get_inner_end,
    $parser->get_outer_end,
  );

  # parse alternative alleles from breakends
  my $parsed_allele;
  my $allele_string;
  if ($so_term =~ /break/ ){
    $parsed_allele = $self->_parse_breakend_alt_alleles($alt, $info);
    $allele_string = $alt;
  }

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
    _line          => $record
  });
  $svf->{allele_string}  = $allele_string if defined $allele_string;
  $svf->{_parsed_allele} = $parsed_allele if defined $parsed_allele;
  $svf->{vep_skip}       = $skip_line     if defined $skip_line;

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
