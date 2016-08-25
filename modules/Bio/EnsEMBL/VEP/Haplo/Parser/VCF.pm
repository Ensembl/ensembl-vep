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

# EnsEMBL module for Bio::EnsEMBL::VEP::Haplo::Parser::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::Parser::VCF - VCF input parser

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::Parser::VCF;

use base qw(Bio::EnsEMBL::VEP::Parser::VCF);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


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

    while(!@$cache) {
      push @$cache, @{$self->create_VariationFeatures()};
      $self->parser->next;
    }
  }

  return @$cache ? shift @$cache : undef;
}

sub create_VariationFeatures {
  my $self = shift;

  my $parser = $self->parser();

  my $gts = $parser->get_samples_genotypes(undef, 1);
  return [] unless scalar keys %$gts;

  return [{
    chr => $parser->get_seqname,
    start => $parser->get_raw_start,
    end => $parser->get_raw_end,
    ids => $parser->get_IDs,
    gts => $gts,
    alleles => $parser->get_reference.','.$parser->get_raw_alternatives,
  }];
}

1;
