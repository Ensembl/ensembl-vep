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

# EnsEMBL module for Bio::EnsEMBL::VEP::Haplo::Parser::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::Parser::VCF - VCF input parser

=head1 SYNOPSIS

my $parser = Bio::EnsEMBL::VEP::Haplo::Parser::VCF->new({
  config            => $config,
  file              => $filename_or_handle,
  valid_chromosomes => $listref_of_chromosome_names,
});

my $vf_hash = $parser->next();

=head1 DESCRIPTION

The Haplo VCF parser class inherits from Bio::EnsEMBL::VEP::Parser::VCF.

It differs from a regular VEP parser in that VariationFeature objects
are not created here, only hash references containing data allowing
them to be lazily created later.

VCF entries must have sample genotypes.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::Parser::VCF;

use base qw(Bio::EnsEMBL::VEP::Parser::VCF);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Variation::Individual;

use Scalar::Util qw(looks_like_number);

=head2 parser

  Example    : $io_parser = $parser->parser;
  Description: Gets the ensembl-io parser object that is used to read data
               from the input file.
  Returntype : Bio::EnsEMBL::IO::Parser::VCF4
  Exceptions : none
  Caller     : all subs
  Status     : Stable

=cut

sub parser {
  my $self = shift;

  if(!exists($self->{parser})) {
    $self->{parser} = $self->SUPER::parser();
    $self->headers();
    $self->samples();
  }

  return $self->{parser};
}


=head2 samples

  Example    : $samples = $parser->samples();
  Description: Gets all sample objects associated with the samples
               represented in the VCF input.
  Returntype : listref of Bio::EnsEMBL::Variation::Sample
  Exceptions : throws if VCF has no sample data
  Caller     : Bio::EnsEMBL::VEP::AnnotationSource classes
  Status     : Stable

=cut

sub samples {
  my $self = shift;

  if(!exists($self->{_samples})) {
    my @sample_ids = @{$self->parser->get_samples};

    throw("ERROR: VCF file contains no sample/individual entries\n") unless @sample_ids;

    $self->{_samples}->{$_} ||= Bio::EnsEMBL::Variation::Sample->new_fast({
      name            => $_,
      display         => 'UNDISPLAYABLE',
      dbID            => --($self->{_sample_id}),
      individual      => Bio::EnsEMBL::Variation::Individual->new_fast({
        name     => $_,
        type_individual => 'outbred',
        dbID     => --($self->{_ind_id}),
      }),
    }) for @sample_ids;
  }

  return $self->{_samples};
}


=head2 next

  Example    : $vf_hash = $parser->next();
  Description: Gets next hash of variant data from input
  Returntype : hashref
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::InputBuffer
  Status     : Stable

=cut

sub next {
  my $self = shift;

  my $vf;
  my $parser = $self->parser;

  while(!$vf) {
    if($self->{_have_read_next}) {
      delete $self->{_have_read_next};
    }
    else {
      $parser->next;
    }
    last unless $parser->{record};
    ($vf) = @{$self->create_VariationFeatures()};
    $self->validate_chr($vf);
  }

  return $vf;
}


=head2 create_VariationFeatures

  Example    : my $vf_hashes = $parser->create_VariationFeatures();
  Description: Gets one or more hashrefs of variant data from input
  Returntype : listref of hashrefs
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub create_VariationFeatures {
  my $self = shift;

  my $parser = $self->parser();

  return [{
    chr => $parser->get_seqname,
    start => $parser->get_raw_start,
    end => $parser->get_raw_end,
    ids => $parser->get_IDs,
    record => \@{$parser->{record}},
    alleles => $parser->get_reference.','.$parser->get_raw_alternatives,
  }];
}

=head2 validate_chr

  Arg 1      : hash
  Example    : my $vf_hashes = $parser->validate_chr($vf);
  Description: Validates chromosome and genomic positions
  Returntype : bool
  Exceptions : none
  Status     : Stable

=cut

sub validate_chr {
  my $self = shift; 
  my $vf   = shift;

  my $vf_chr = $vf->{chr}; 

  # sanity checks
  unless(looks_like_number($vf->{start}) && looks_like_number($vf->{end})) {
    $self->warning_msg("WARNING: Start ".$vf->{start}." or end ".$vf->{end}." coordinate invalid on line ".$self->line_number);
    return 0;
  } 

  unless($self->_have_chr($vf)){
    $self->warning_msg("Chromosome ".$vf->{chr}." not found in annotation sources on line ".$self->line_number);
    return 0;
  } 

  return 1;
}

1;
