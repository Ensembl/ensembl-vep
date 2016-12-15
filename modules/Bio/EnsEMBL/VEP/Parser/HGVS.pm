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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::HGVS
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::HGVS - HGVS list input parser

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Parser::HGVS;

use base qw(Bio::EnsEMBL::VEP::Parser);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::ListBasedParser;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # requires db connection
  throw("ERROR: Cannot use HGVS format in offline mode") if $self->param('offline');

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

  my $hgvs = $parser->get_value;

  # remove whitespace
  $hgvs =~ s/\s+//g;

  my $core_group = $self->param('core_type');

  $DB::single = 1;
  
  my $vfa = $self->get_adaptor('variation', 'VariationFeature');
  my $sa  = $self->get_adaptor($core_group, 'Slice');
  my $ta  = $self->get_adaptor($core_group, 'Transcript');

  my $vf;

  # not all hgvs notations are supported yet, so we have to wrap it in an eval
  eval { $vf = $vfa->fetch_by_hgvs_notation($hgvs, $sa, $ta) };

  if(!defined($vf) || (defined $@ && length($@) > 1)) {
    $self->warning_msg("WARNING: Unable to parse HGVS notation \'$hgvs\'\n$@");
    return [];
  }

  # transfer to whole chromosome slice
  $vf = $vf->transfer($vf->slice->seq_region_Slice);

  # name it after the HGVS
  $vf->{variation_name} = $hgvs;

  # add chr attrib
  $vf->{chr} = $vf->slice->seq_region_name;

  $vf->{_line} = [$hgvs];

  return $self->post_process_vfs([$vf]);
}

1;