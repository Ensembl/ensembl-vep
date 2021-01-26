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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::VEP_input
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::VEP_input - VEP format input parser

=head1 SYNOPSIS

my $parser = Bio::EnsEMBL::VEP::Parser::VEP_input->new({
  config => $config,
  file   => 'variants.txt',
});

my $vf = $parser->next();

=head1 DESCRIPTION

VEP format input parser.

A simple tab-delimited format to describe sequence or simple structural
variants.

http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default

=head1 METHODS

=cut


use strict;
use warnings;
no warnings 'recursion';

package Bio::EnsEMBL::VEP::Parser::VEP_input;

use base qw(Bio::EnsEMBL::VEP::Parser);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::Parser::VEP_input;


=head2 parser

  Example    : $io_parser = $parser->parser();
  Description: Get ensembl-io parser object used to read data from input.
  Returntype : Bio::EnsEMBL::IO::Parser::VEP_input
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;

  if(!exists($self->{parser})) {
    $self->{parser} = Bio::EnsEMBL::IO::Parser::VEP_input->open($self->file);
    $self->{parser}->{delimiter} = $self->delimiter;
  }

  return $self->{parser};
}


=head2 create_VariationFeatures

  Example    : $vfs = $parser->create_VariationFeatures();
  Description: Create a VariationFeature object from the current line
               of input. 
  Returntype : arrayref of Bio::EnsEMBL::BaseVariationFeature
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub create_VariationFeatures {
  my $self = shift;

  my $parser = $self->parser;
  $parser->next();

  $self->skip_empty_lines();

  return [] unless $parser->{record};

  $self->line_number($self->line_number + 1);

  my ($chr, $start, $end, $allele_string, $strand, $var_name) = map {$_ =~ s/\s+$// if defined($_); $_} (
    $parser->get_seqname,
    $parser->get_start,
    $parser->get_end,
    $parser->get_allele,
    $parser->get_strand,
    $parser->get_id
  );

  # check strand
  $strand = ($strand || '') =~ /\-/ ? -1 : 1;

  my $vf;

  # sv?
  if($allele_string !~ /\//) {
    my $so_term;

    # convert to SO term
    my %terms = (
      INS  => 'insertion',
      DEL  => 'deletion',
      TDUP => 'tandem_duplication',
      DUP  => 'duplication'
    );

    $so_term = defined $terms{$allele_string} ? $terms{$allele_string} : $allele_string;

    $vf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
      start          => $start,
      end            => $end,
      strand         => $strand,
      adaptor        => $self->get_adaptor('variation', 'StructuralVariationFeature'),
      variation_name => $var_name,
      chr            => $chr,
      class_SO_term  => $so_term,
    });
  }

  # normal vf
  else {
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
      start          => $start,
      end            => $end,
      allele_string  => $allele_string,
      strand         => $strand,
      map_weight     => 1,
      adaptor        => $self->get_adaptor('variation', 'VariationFeature'),
      variation_name => $var_name,
      chr            => $chr,
    });
  }

  $vf->{_line} = $parser->{record};

  return $self->post_process_vfs([$vf]);
}

1;
