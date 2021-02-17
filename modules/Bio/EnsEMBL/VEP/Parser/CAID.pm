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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::CAID
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::CAID - CAid list input parser

=head1 SYNOPSIS

my $parser = Bio::EnsEMBL::VEP::Parser::CAID->new({
  config => $config,
  file   => 'variant_caids.txt',
});

my $vf = $parser->next();

=head1 DESCRIPTION

CAid - ClinGen Canonical Allele identifier format parser.

See:  ClinGen Allele Registry: https://doi.org/10.1002/humu.23637

CAids are looked up in the Ensembl homo_sapiens Variation database.

CAids stored in the database are downloaded from http://reg.genome.network

One CAid may return multiple VariationFeature objects as several variants may
have equivalent alleles or be co-located. Therefore as for similar parsers
a method create_VariationFeatures returns an array ref of VariationFeatures.


=head1 METHODS

=cut


use strict;
use warnings;
no warnings 'recursion';

package Bio::EnsEMBL::VEP::Parser::CAID;

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
  Example    : $parser = Bio::EnsEMBL::VEP::Parser::CAID->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::Parser::CAID object.
  Returntype : Bio::EnsEMBL::VEP::Parser::CAID
  Exceptions : throws if offline mode (--offline) enabled
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  # requires db connection
  throw("ERROR: Cannot use CAid format in offline mode") if $self->param('offline');

  $self->{adaptor} = $self->get_adaptor('variation', 'VariationFeature');

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
  Description: Create one or more VariationFeature objects from the current line
               of input; multiple may be returned if a CAid maps to multiple
               variant alleles
  Returntype : arrayref of Bio::EnsEMBL::VariationFeature
  Exceptions : warns if CAid not found in database
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

  my $caid = $parser->get_value;

  # remove whitespace
  $caid =~ s/\s+//g;

  my $ad = $self->{adaptor};

  # tell adaptor to fetch failed variants
  # but store state to restore afterwards
  my $prev = $ad->db->include_failed_variations;
  $ad->db->include_failed_variations(1);
  my @vfs = @{$ad->fetch_all_by_caid($caid)};

  unless (@vfs) {
    $self->warning_msg("WARNING: No variants found for CAid \'$caid\'");
    return $self->create_VariationFeatures();
  }

  for (@vfs) {
    delete $_->{dbID};
    delete $_->{overlap_consequences};
    $_->{chr} = $_->seq_region_name;
    $_->{variation_name} = $caid;
    $_->{_line} = [$caid];
  }

  # restore state
  $ad->db->include_failed_variations($prev);
  return $self->post_process_vfs(\@vfs);
}

1;
