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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::VCF - VCF input parser

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Parser::VCF;

use parent qw(Bio::EnsEMBL::VEP::Parser);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::Parser::VCF4;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(allow_non_variant gp individual process_ref_homs phased)]);

  return $self;
}

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::VCF4->open($self->file);
}

sub next {
  my $self = shift;

  my $cache = $self->{_vf_cache} ||= [];

  if(!scalar @$cache) {
    $self->parser->next;
    push @$cache, @{$self->create_VariationFeatures()};

    $self->line_number($self->line_number + 1) if scalar @$cache;
  }

  return shift @$cache;
}

sub create_VariationFeatures {
  my $self = shift;

  my $parser = $self->parser();
  my $record = $parser->{record};

  return [] unless $record;

  # get the data we need to decide if this is an SV
  my ($alts, $info) = (
    $parser->get_alternatives,
    $parser->get_info,
  );

  if($info->{SVTYPE} || join(",", @$alts) =~ /[<\[][^\*]+[>\]]/) {
    return $self->create_StructuralVariationFeatures();
  }

  # get the rest of the relevant data
  my ($chr, $start, $end, $ref, $ids) = (
    $parser->get_seqname,
    $parser->get_raw_start,
    $parser->get_raw_end,
    $parser->get_reference,
    $parser->get_IDs,
  );

  # non-variant
  my $non_variant = 0;

  if($alts->[0] eq '.') {
    if($self->{allow_non_variant}) {
      $non_variant = 1;
    }
    else {
      return [];
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
      $self->warning_msg("No GP flag found in INFO column on line ".$self->line_number);
      return [];
    }
  }

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
  });

  # flag as non-variant
  $vf->{non_variant} = 1 if $non_variant;

  # individual data?
  return $self->create_individual_VariationFeatures($vf) if $self->{individual};

  # normal return
  return [$vf];
}

sub create_StructuralVariationFeatures {
  my $self = shift;

  my $parser = $self->parser();
  my $record = $parser->{record};

  # get relevant data
  my ($chr, $start, $end, $alts, $info, $ids) = (
    $parser->get_seqname,
    $parser->get_start,
    $parser->get_end,
    $parser->get_alternatives,
    $parser->get_info,
    $parser->get_IDs,
  );

  my $alt = join(",", @$alts);

  # work out the end coord
  if(defined($info->{END})) {
    $end = $info->{END};
  }
  elsif(defined($info->{SVLEN})) {
    $end = $start + abs($info->{SVLEN}) - 1;
  }

  # check for imprecise breakpoints
  my ($min_start, $max_start, $min_end, $max_end) = (
    $parser->get_outer_start,
    $parser->get_inner_start,
    $parser->get_inner_end,
    $parser->get_outer_end,
  );

  # get type
  my $type;

  if($alt =~ /\<|\[|\]|\>/) {
    $type = $alt;
    $type =~ s/\<|\>//g;
    $type =~ s/\:.+//g;

    if($start >= $end && $type =~ /del/i) {
      my $line = join("\t", @$record);
      $self->warning_msg("WARNING: VCF line on line ".$self->line_number." looks incomplete, skipping:\n$line\n");
      return [];
    }

  }
  else {
    $type = $info->{SVTYPE};
  }

  my $so_term;

  if(defined($type)) {
    # convert to SO term
    my %terms = (
      INS  => 'insertion',
      DEL  => 'deletion',
      TDUP => 'tandem_duplication',
      DUP  => 'duplication'
    );

    $so_term = defined $terms{$type} ? $terms{$type} : $type;
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
    chr            => $chr,
    class_SO_term  => $so_term,
  });

  return [$svf];
}

sub create_individual_VariationFeatures {
  my $self = shift;
  my $vf = shift;

  my $parser = $self->parser();
  my $record = $parser->{record};

  my @alleles = split '\/', $vf->{allele_string};
  my $ref = $alleles[0];

  my @return;

  # get genotypes from parser
  my $include = lc($self->{individual}->[0]) eq 'all' ? $parser->get_samples : $self->{individual};
  my $ind_gts = $parser->get_samples_genotypes($include, 1 - ($self->{allow_non_variant} || 0));

  foreach my $ind(@$include) {

    # get alleles present in this individual
    my $gt = $ind_gts->{$ind};
    my @bits = split /\||\/|\\/, $gt;
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

1;
