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
use Bio::EnsEMBL::VEP::Utils qw(trim_sequences);

our $HASH_TREE_SIZE = 1e4;
our $CAN_USE_INTERVAL_TREE;

BEGIN {
  if (eval q{ require Set::IntervalTree; 1 }) {
    $CAN_USE_INTERVAL_TREE = 1;
  }
}

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(buffer_size minimal)]);

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
    if($prev_chr && $vf->{chr} ne $prev_chr && $vf->{chr} !~ /LRG/ && $prev_chr !~ /LRG/) {
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
      if($prev_chr && $vf->{chr} ne $prev_chr && $vf->{chr} !~ /LRG/ && $prev_chr !~ /LRG/) {

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

  $self->split_variants() if $self->{minimal};

  return $buffer;
}

sub get_overlapping_vfs {
  my $self = shift;
  my $start = shift;
  my $end = shift;

  ($start, $end) = ($end, $start) if $start > $end;

  if(my $tree = $self->interval_tree) {
    return [
      grep { overlap($_->{start}, $_->{end}, $start, $end) }
      @{$tree->fetch($start - 1, $end)}
    ];
  }
  else {
    my $hash_tree = $self->hash_tree;

    return [
      grep { overlap($_->{start}, $_->{end}, $start, $end) }
      values %{{
        map {$_->{_hash_tree_id} => $_}                  # use _hash_tree_id to uniquify
        map {@{$hash_tree->{$_} || []}}                  # tree might be empty
        (
          int($start / $HASH_TREE_SIZE)      # start of range
          ..
          (int($end / $HASH_TREE_SIZE) + 1)  # end of range
        )
      }}
    ];
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

# this is a pure-perl interval tree of sorts
# variants are added to "branches" (keys) of a hash
# where the key is a bin obtained by bin = int(coord / $HASH_TREE_SIZE)
# variants can span multiple bins and so can be added to more than one branch
sub hash_tree {
  my $self = shift;

  if(!exists($self->{temp}->{hash_tree})) {

    my $hash_tree = {};
    my $hash_tree_id = 1;

    foreach my $vf(@{$self->buffer}) {
      my ($vf_s, $vf_e) = ($vf->{start}, $vf->{end});
      ($vf_s, $vf_e) = ($vf_e, $vf_s) if $vf_s > $vf_e;
      my ($h_s, $h_e) = map {int($_ / $HASH_TREE_SIZE)} ($vf_s, $vf_e);

      # add this VF to the bin for each branch that it overlaps
      for(my $s = $h_s; $s <= $h_e; $s++) {
        push @{$hash_tree->{$s}}, $vf;
      }

      # log a unique ID for each VF in the tree
      # allows us to unique sort when fetching across multiple branches
      $vf->{_hash_tree_id} = $hash_tree_id++;
    }

    $self->{temp}->{hash_tree} = $hash_tree;
  }

  return $self->{temp}->{hash_tree};
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
  $self->stats->log_lines_read($self->parser->line_number) if $self->parser;
  
  foreach my $vf(@{$self->buffer}) {
    $vf->{slice} ||= $self->get_slice($vf->{chr});
    $vf->_finish_annotation;
  }
}

sub reset_buffer {
  $_[0]->{temp} = {};
}

sub reset_pre_buffer {
  $_[0]->{pre_buffer} = [];
}

sub buffer {
  my $self = shift;
  $self->{temp}->{buffer} = shift if @_;
  return $self->{temp}->{buffer} ||= [];
}

sub pre_buffer {
  return $_[0]->{pre_buffer} ||= [];
}

sub rejoin_required {
  my $self = shift;
  $self->{temp}->{rejoin_required} = shift if @_;
  return $self->{temp}->{rejoin_required} ||= 0;
}

sub split_variants {
  my $self = shift;
  my $listref = $self->buffer;

  my $original_vf_count = scalar @$listref;

  # split and link multi-allele VFs
  my @split_list;

  foreach my $original_vf(@$listref)  {
    if($original_vf->{allele_string} =~ /.+\/.+\/.+/) {

      my @alleles = split('/', $original_vf->{allele_string});
      my $original_ref = shift @alleles;
      my $first;

      my @tmp;
      my $changed = 0;
      my $base_allele_number = 1;

      foreach my $alt(@alleles) {

        my $ref   = $original_ref;
        my $start = $original_vf->{start};
        my $end   = $original_vf->{end};
        my $this_changed = 0;

        ($ref, $alt, $start, $end, $this_changed) = @{trim_sequences($ref, $alt, $start, $end)};
        $ref ||= '-';
        $alt ||= '-';
        $changed += $this_changed;

        # create a copy
        my $new_vf;
        %$new_vf = %{$original_vf};
        bless $new_vf, ref($original_vf);

        # give it a new allele string and coords
        $new_vf->allele_string($ref.'/'.$alt);
        $new_vf->{start} = $start;
        $new_vf->{end} = $end;
        $new_vf->{alt_allele} = $alt;

        # $new_vf->{variation_name} = 'merge_'.$alt;

        # not the first one ($first already exists)
        if($first) {
          $new_vf->{merge_with} = $first;
          $new_vf->{_base_allele_number} = $base_allele_number++;
        }

        # this is the first one
        else {
          $first = $new_vf;

          # $new_vf->{variation_name} = 'first_'.$alt;

          # store the original allele string and coords
          $first->{original_allele_string} = $original_vf->{allele_string};
          $first->{original_start}         = $original_vf->{start};
          $first->{original_end}           = $original_vf->{end};
          $first->{minimised}              = 1
        }

        push @tmp, $new_vf;
      }

      if($changed) {
        push @split_list, @tmp;
      }
      else {
        push @split_list, $original_vf;
      }
    }
    else {
      push @split_list, $original_vf;
    }
  }

  $self->buffer(\@split_list);
  $self->rejoin_required(scalar @split_list != $original_vf_count);
}

1;