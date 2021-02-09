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

# EnsEMBL module for Bio::EnsEMBL::VEP::Haplo::InputBuffer
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::InputBuffer - class representing a buffer of VariationFeatures to be processed

=head1 SYNOPSIS

my $ib = Bio::EnsEMBL::VEP::Haplo::InputBuffer->new({
  config => $config,
  parser => $parser,
  transcript_tree => $tree
});

$ib->next();

=head1 DESCRIPTION

The Haplo InputBuffer class inherits from the Bio::EnsEMBL::VEP::InputBuffer,
and offers similar functionality. It is intended to contain a buffer of
VariationFeature objects as read from a Parser.

It differs from its parent class in how the buffer is filled. On reading the
first variant from the parser, the transcript tree is used to find any
overlapping transcripts with this variant; the buffer is then filled by reading
in variants until the coordinates no longer overlap the transcripts'.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::InputBuffer;

use base qw(Bio::EnsEMBL::VEP::InputBuffer);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 new

  Arg 1      : hashref $args
               {
                 config          => Bio::EnsEMBL::VEP::Config,
                 parser          => Bio::EnsEMBL::VEP::Haplo::Parser::VCF,
                 transcript_tree => Bio::EnsEMBL::VEP::TranscriptTree,
               }
  Example    : $ib = Bio::EnsEMBL::VEP::Haplo::InputBuffer->new($args);
  Description: Creates an InputBuffer object.
  Returntype : Bio::EnsEMBL::VEP::Haplo::InputBuffer
  Exceptions : throws if no transcript tree or wrong class
  Caller     : Bio::EnsEMBL::VEP::Haplo::Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  my $tree = $hashref->{transcript_tree};
  assert_ref($tree, 'Bio::EnsEMBL::VEP::TranscriptTree');
  $self->transcript_tree($tree);

  return $self;
}


=head2 next

  Example    : $ib->next()
  Description: Resets buffer and fills with variants overlapping next found transcript(s)
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature objects
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::Haplo::Runner
  Status     : Stable

=cut

sub next {
  my $self = shift;

  $self->reset_buffer();

  my $pre_buffer = $self->pre_buffer();
  my $buffer = $self->buffer();

  if(my $parser = $self->parser) {
    my $max = 0;
    my $current_chr;
    my $vf;

    # Setup new buffer max
    while(!$max && ($vf = @$pre_buffer ? shift @$pre_buffer : $parser->next)) {
      $max = $self->get_max_from_tree($vf->{chr}, $vf->{start}, $vf->{end});
      $current_chr = $vf->{chr};
    }

    return $buffer unless $max;

    # Continue to add VF to buffer until VF start > max or VF on different chromosome
    while($vf && $vf->{start} <= $max && $vf->{chr} eq $current_chr) {
      push @$buffer, $vf;
      $vf = $parser->next;
    }

    push @$pre_buffer, $vf if $vf;
  }

  return $buffer;
}


=head2 get_max_from_tree
  
  Arg 1      : string $chromosome
  Arg 2      : int $start
  Arg 3      : int $end
  Example    : my $max = $ib->get_max_from_tree($chr, $start, $end)
  Description: Given a location, find the end of the block of transcripts
               that overlap. Essentially finds the beginning of the next
               intergenic region. Returns 0 if no overlapping transcripts.
  Returntype : int
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub get_max_from_tree {
  my ($self, $c, $s, $e) = @_;

  my $tree = $self->transcript_tree;

  my $prev_max = 0;

  while(1) {
    my @regions = @{$tree->fetch($c, $s - 1, $e)};
    last unless @regions;

    foreach my $r(@regions) {
      ($s, $e) = @$r if $r->[1] > $e;
    }

    last if $prev_max == $e;
    $prev_max = $e;
  }

  return $prev_max;
}


=head2 transcript_tree
  
  Arg 1      : (optional) Bio::EnsEMBL::VEP::TranscriptTree $tree
  Example    : my $tree = $ib->transcript_tree
  Description: Getter/setter for the transcript tree object associated with this
               input buffer.
  Returntype : Bio::EnsEMBL::VEP::TranscriptTree
  Exceptions : none
  Caller     : get_max_from_tree()
  Status     : Stable

=cut

sub transcript_tree {
  my $self = shift;
  $self->{transcript_tree} = shift if @_;
  return $self->{transcript_tree};
}

1;

