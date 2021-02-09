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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser::SPDI
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser::SPDI - SPDI list input parser

=head1 SYNOPSIS

my $parser = Bio::EnsEMBL::VEP::Parser::SPDI->new({
  config => $config,
  file   => 'spdi.txt',
});

my $vf = $parser->next();

=head1 DESCRIPTION

SPDI format parser.

See https://www.ncbi.nlm.nih.gov/variation/notation/ for spec.

Variants can be genomic.

=head1 METHODS

=cut


use strict;
use warnings;
no warnings 'recursion';

package Bio::EnsEMBL::VEP::Parser::SPDI;

use base qw(Bio::EnsEMBL::VEP::Parser);

use Bio::EnsEMBL::VEP::Utils qw(merge_arrays);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::IO::ListBasedParser;

=head2 new

  Arg 1      : hashref $args
               {
                 config    => Bio::EnsEMBL::VEP::Config,
                 file      => string or filehandle,
               }
  Example    : $parser = Bio::EnsEMBL::VEP::Parser::SPDI->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::Parser::SPDI object.
  Returntype : Bio::EnsEMBL::VEP::Parser::SPDI
  Exceptions : throws if offline mode (--offline) enabled
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  # requires db connection
  throw("ERROR: Cannot use SPDI format in offline mode") if $self->param('offline');

  $self->add_shortcuts(['ambiguous_spdi']);

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
  Exceptions : warns if unable to parse SPDI string
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

  my $spdi = $parser->get_value;

  # remove whitespace
  $spdi =~ s/\s+//g;

  my $param_core_group = $self->param('core_type');

  my @core_groups = sort {($b eq $param_core_group) cmp ($a eq $param_core_group)} qw(core otherfeatures);

  my $vfa = $self->get_adaptor('variation', 'VariationFeature');
  my $vfs = [];
  my @errors;

  foreach my $core_group(@core_groups) {
    my $sa  = $self->get_adaptor($core_group, 'Slice');
    my $ta  = $self->get_adaptor($core_group, 'Transcript');

    eval {
        push @$vfs, $vfa->fetch_by_spdi_notation(
          -spdi               => $spdi,
          -slice_adaptor      => $sa,
          -transcript_adaptor => $ta,
          -replace_ref        => $self->{lookup_ref} || 0,
        );
    };

    # only log unique errors
    merge_arrays(\@errors, [$@]) if $@ && length($@) > 1;

    last if @$vfs;
  }

  unless(@$vfs || scalar(@errors) == 0) {
    my %known_messages_hash = ('MSG: Region requested must be smaller than 5kb' => 0);

    my @grep_names = grep(/^MSG:/, split(/\n/, $errors[0]));
    my @error_message = exists( $known_messages_hash{$grep_names[0]}) ? @grep_names : @errors;

    $self->warning_msg("WARNING: Unable to parse SPDI notation \'$spdi\'\n".join("\n", @error_message));
    return $self->create_VariationFeatures;
  }

  foreach my $vf(@$vfs) {

    # transfer to whole chromosome slice
    $vf = $vf->transfer($vf->slice->seq_region_Slice);

    # name it after the SPDI
    $vf->{variation_name} = $spdi;

    # add chr attrib
    $vf->{chr} = $vf->slice->seq_region_name;

    $vf->{_line} = [$spdi];
  }

  # post_process_vfs will do lookup_ref again, we've already done it
  my $prev = $self->{lookup_ref} || 0;
  $self->{lookup_ref} = 0;
  my $return = $self->post_process_vfs($vfs);
  $self->{lookup_ref} = $prev;

  return $return;
}

1;
