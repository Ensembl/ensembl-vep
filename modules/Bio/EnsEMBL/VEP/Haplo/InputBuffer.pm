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

# EnsEMBL module for Bio::EnsEMBL::VEP::Haplo::InputBuffer
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::InputBuffer - class representing a buffer of VariationFeatures to be processed

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::InputBuffer;

use base qw(Bio::EnsEMBL::VEP::InputBuffer);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  my $tree = $hashref->{transcript_tree};
  assert_ref($tree, 'Bio::EnsEMBL::VEP::Haplo::TranscriptTree');
  $self->transcript_tree($tree);

  return $self;
}

sub next {
  my $self = shift;

  $self->reset_buffer();

  my $pre_buffer = $self->pre_buffer();
  my $buffer = $self->buffer();

  if(my $parser = $self->parser) {

    my $max = 0;
    my $vf;

    while(!$max && ($vf = $parser->next)) {
      $max = $self->get_max_from_tree($vf->{chr}, $vf->{start}, $vf->{end});
    }

    return $buffer unless $max;

    while($vf && $vf->{start} <= $max) {
      push @$buffer, $vf;
      $vf = $parser->next;
    }
  }

  return $buffer;
}

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

sub transcript_tree {
  my $self = shift;
  $self->{transcript_tree} = shift if @_;
  return $self->{transcript_tree};
}

1;