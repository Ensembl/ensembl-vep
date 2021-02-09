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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::Region
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::Region - region format parser

=head1 SYNOPSIS

my $parser = Bio::EnsEMBL::VEP::Parser::Region->new({
  config => $config,
  file   => 'regions.txt',
});

my $vf = $parser->next();

=head1 DESCRIPTION

Region format parser. This is a format included to support
REST API lookups that use the format e.g. REST URL:

http://rest.ensembl.org/vep/human/region/9:22125503-22125502:1/C

Only the final two "/"-separated components form the format.

[chr]:[start]-[end]:[strand]/[allele]

chr, start, end and allele are mandatory fields.

chr, start and end denote the start and end of the region of DNA
to be replaced by the given allele. To represent an insertion
between pos X and pos X+1, use [chr]:[X+1]-[X]:[strand]/[allele]

strand is optional, if absent it defaults to a value of 1 (forward) strand.
Valid values are "1", "+1", "-1".

allele may be any DNA string; "-" denotes a deletion.

"DUP" and "DEL" may be used as alleles to denote larger structural
variants.

=head1 METHODS

=cut


use strict;
use warnings;
no warnings 'recursion';

package Bio::EnsEMBL::VEP::Parser::Region;

use base qw(Bio::EnsEMBL::VEP::Parser);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::ListBasedParser;


=head2 new

  Arg 1      : hashref $args
               {
                 config    => Bio::EnsEMBL::VEP::Config,
                 file      => string or filehandle,
               }
  Example    : $parser = Bio::EnsEMBL::VEP::Parser::Region->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::Parser::Region object.
  Returntype : Bio::EnsEMBL::VEP::Parser::Region
  Exceptions : throws if no database or FASTA available for reference allele lookup
               throws if --check_ref param set
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # requires db connection
  throw("ERROR: Cannot use Region format in offline mode without a FASTA file") if $self->param('offline') && !$self->param('fasta');

  # can't be used with --check_ref as no ref given
  throw("ERROR: Region format is not compatible with --check_ref") if $self->param('check_ref');

  # config so lookups get done
  $self->param('lookup_ref', 1);
  $self->{lookup_ref} = 1;

  return $self;
}


=head2 parser

  Example    : $io_parser = $parser->parser();
  Description: Get ensembl-io parser object used to read data from input.
  Returntype : Bio::EnsEMBL::IO::ListBasedParser
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::ListBasedParser->open($self->file);
}


=head2 create_VariationFeatures

  Example    : $vfs = $parser->create_VariationFeatures();
  Description: Create a VariationFeature object from the current line
               of input. 
  Returntype : arrayref of Bio::EnsEMBL::VariationFeature
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

  my $region = $parser->get_value();

  return [] unless $region =~ /^([^\:]+)\:(\d+)\-(\d+)(\:[\-\+]?1)?[\/\:](ins|dup|del|[ACGTN-]+)$/i;
  my ($chr, $start, $end) = ($1, $2, $3);

  my ($strand, $allele);
  if($5) {
    ($strand, $allele) = ($4, uc($5));
  }
  else {
    $allele = uc($4);
  }

  # check strand
  $strand = ($strand || '') =~ /\-/ ? -1 : 1;

  my $vf;
  
  # sv?
  if($allele =~ /^(INS|DEL|DUP)$/) {
    my $so_term;

    # convert to SO term
    my %terms = (
      INS  => 'insertion',
      DEL  => 'deletion',
      TDUP => 'tandem_duplication',
      DUP  => 'duplication'
    );

    $so_term = defined $terms{$allele} ? $terms{$allele} : $allele;

    $vf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
      start          => $start,
      end            => $end,
      strand         => $strand,
      adaptor        => $self->get_adaptor('variation', 'StructuralVariationFeature'),
      variation_name => $region,
      chr            => $chr,
      class_SO_term  => $so_term,
    });
  }

  # normal vf
  else {
    my $ref = ('N' x (($end - $start) + 1)) || '-';

    $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
      start          => $start,
      end            => $end,
      allele_string  => $ref.'/'.$allele,
      strand         => $strand,
      map_weight     => 1,
      adaptor        => $self->get_adaptor('variation', 'VariationFeature'),
      variation_name => $region,
      chr            => $chr,
    });
  }

  $vf->{_line} = [$region];

  return $self->post_process_vfs([$vf]);
}

1;