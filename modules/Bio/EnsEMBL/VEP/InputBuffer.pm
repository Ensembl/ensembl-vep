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

# EnsEMBL module for Bio::EnsEMBL::VEP::InputBuffer
#
#

=head1 NAME

Bio::EnsEMBL::VEP::InputBuffer - class representing a buffer of VariationFeatures to be processed by VEP

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::InputBuffer;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

our $CAN_USE_INTERVAL_TREE;

BEGIN {
  if (eval { require Set::IntervalTree; 1 }) {
    $CAN_USE_INTERVAL_TREE = 1;
  }
  else {
    $CAN_USE_INTERVAL_TREE = 0;
  }
}

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(buffer_size)]);

  my $hashref = $_[0];
  if($hashref) {
    $self->parser($hashref->{parser}) if $hashref->{parser};

    if($hashref->{variation_features}) {
      my $buffer = $self->pre_buffer;
      push @$buffer, @{$hashref->{variation_features}};
    }
  }

  return $self;
}

sub parser {
  my $self = shift;

  if(@_) {
    my $parser = shift;
    assert_ref($parser, 'Bio::EnsEMBL::VEP::Parser');
    $self->{parser} = $parser;
  }

  return $self->{parser};
}

sub next {
  my $self = shift;

  $self->reset_buffer();

  my $pre_buffer = $self->pre_buffer();
  my $buffer = $self->buffer();

  my $buffer_size = $self->{buffer_size};

  # We do a bit of trickery here to help memory usage later.
  # Basically we don't want the buffer to contain variants
  # from more than one chromosome, so we "shortcut" out
  # filling the buffer before it hits $buffer_size.
  my $prev_chr;

  while(@$pre_buffer && @$buffer < $buffer_size) {
    my $vf = $pre_buffer->[0];
    
    # new chromosome
    if($prev_chr && $vf->{chr} ne $prev_chr) {
      return $buffer;
    }

    # same chromosome
    else {
      push @$buffer, shift @$pre_buffer;
      $prev_chr = $vf->{chr};
    }
  }

  if(my $parser = $self->parser) {
    while(@$buffer < $buffer_size && (my $vf = $parser->next)) {

      # new chromosome
      if($prev_chr && $vf->{chr} ne $prev_chr) {

        # we can't push the VF back onto the parser, so add it to $pre_buffer
        # and it will get picked up on the following next() call
        push @$pre_buffer, $vf;
        return $buffer;
      }

      # same chromosome
      else {
        push @$buffer, $vf;
        $prev_chr = $vf->{chr};
      }
    }
  }

  return $buffer;
}

sub get_overlapping_vfs {
  my $self = shift;
  my $start = shift;
  my $end = shift;

  ($start, $end) = ($end, $start) if $start > $end;

  if(my $tree = $self->interval_tree) {
    return $tree->fetch($start - 1, $end);
  }
  else {
    return [grep {overlap($_->{start}, $_->{end}, $start, $end)} @{$self->buffer}];
  }
}

sub interval_tree {
  my $self = shift;

  if(!exists($self->{temp}->{interval_tree})) {

    return $self->{temp}->{interval_tree} = undef unless $CAN_USE_INTERVAL_TREE;

    my $tree = Set::IntervalTree->new();

    foreach my $vf(@{$self->buffer}) {
      my ($s, $e) = ($vf->{start}, $vf->{end});
      ($s, $e) = ($e, $s) if $s > $e;
      $tree->insert($vf, $s - 1, $e);
    }

    $self->{temp}->{interval_tree} = $tree;
  }

  return $self->{temp}->{interval_tree};
}

sub min_max {
  my $self = shift;

  if(!exists($self->{temp}->{min_max})) {
    return $self->{temp}->{min_max} = shift if @_;

    my ($min, $max) = (1e10, 0);

    foreach my $vf(@{$self->buffer}) {
      my ($vf_s, $vf_e) = ($vf->{start}, $vf->{end});

      if($vf_s > $vf_e) {
        $min = $vf_e if $vf_e < $min;
        $max = $vf_s if $vf_s > $max;
      }
      else {
        $min = $vf_s if $vf_s < $min;
        $max = $vf_e if $vf_e > $max;
      }
    }

    $self->{temp}->{min_max} = [$min, $max];
  }

  return $self->{temp}->{min_max};
}

sub finish_annotation {
  my $self = shift;
  $_->_finish_annotation for @{$self->buffer};
}

sub reset_buffer {
  $_[0]->{temp} = {};
}

sub buffer {
  return $_[0]->{temp}->{buffer} ||= [];
}

sub pre_buffer {
  return $_[0]->{pre_buffer} ||= [];
}

1;