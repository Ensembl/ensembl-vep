=head1 LICENSE

Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::ID
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::ID - ID list input parser

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Parser::ID;

use base qw(Bio::EnsEMBL::VEP::Parser);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::ListBasedParser;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # requires db connection
  throw("ERROR: Cannot use ID format in offline mode") if $self->param('offline');

  $self->{adaptor} = $self->get_adaptor('variation', 'Variation');

  return $self;
}

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::ListBasedParser->open($self->file);
}

sub create_VariationFeatures {
  my $self = shift;

  my $parser = $self->parser;
  $parser->next();

  $self->skip_empty_lines();

  return [] unless $parser->{record};

  $self->line_number($self->line_number + 1);

  my $id = $parser->get_value;

  # remove whitespace
  $id =~ s/\s+//g;

  my $ad = $self->{adaptor};

  # tell adaptor to fetch failed variants
  # but store state to restore afterwards
  my $prev = $ad->db->include_failed_variations;
  $ad->db->include_failed_variations(1);

  my $v_obj = $ad->fetch_by_name($id);

  unless($v_obj) {
    $self->warning_msg("WARNING: No variant found with ID \'$id\'");
    return $self->create_VariationFeatures();
  }

  my @vfs = @{$v_obj->get_all_VariationFeatures};

  unless(@vfs) {
    $self->warning_msg("WARNING: No mappings found for variant \'$id\'");
    return $self->create_VariationFeatures();
  }

  for(@vfs) {
    delete $_->{dbID};
    delete $_->{overlap_consequences};
    $_->{chr} = $_->seq_region_name;
    $_->{variation_name} = $id;
    $_->{_line} = [$id];
  }

  # restore state
  $ad->db->include_failed_variations($prev);

  return $self->post_process_vfs(\@vfs);
}

1;