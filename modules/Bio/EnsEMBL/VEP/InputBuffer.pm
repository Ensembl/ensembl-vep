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

# EnsEMBL module for Bio::EnsEMBL::VEP::InputBuffer
#
#

=head1 NAME

Bio::EnsEMBL::VEP::InputBuffer - class representing a buffer of VariationFeatures to be processed by VEP

=head1 SYNOPSIS

my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({
  config => $config,
  parser => $parser 
});

my $vfs = $ib->next();

=head1 DESCRIPTION

The InputBuffer class represents a chunk or buffer of VariationFeature
objects as created by the attached Parser object. Calling next() fills
up the buffer to the size set by the parameter buffer_size.

Also contains methods for retrieving VFs in the current buffer that
overlap a given coordinate range - this is done using either
Set::IntervalTree (preferred if installed) or a pure perl hash-based
lookup.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::InputBuffer;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(trim_sequences);


our $HASH_TREE_SIZE = 1e4;
our $CAN_USE_INTERVAL_TREE;

BEGIN {
  if (eval q{ require Set::IntervalTree; 1 }) {
    $CAN_USE_INTERVAL_TREE = 1;
  }
}


=head2 new

  Arg 1      : hashref $args
               {
                 config             => Bio::EnsEMBL::VEP::Config,
                 parser             => Bio::EnsEMBL::VEP::Parser,
                 variation_features => listref of Bio::EnsEMBL::Variation::VariationFeature (optional, added to pre-buffer),
               }
  Example    : $ib = Bio::EnsEMBL::VEP::InputBuffer->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::InputBuffer object.
               If the variation_features arg is supplied, the pre-buffer
               is filled with those VariationFeatures.
  Returntype : Bio::EnsEMBL::VEP::InputBuffer
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

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
    # Possibility to change the default value of the "max_not_ordered_variants"
    # Could also implement a way to change it from the Runner object by fetching 
    # a "$self->param('max_not_ordered_variants')" - for future developments
    $self->{max_not_ordered_variants} = $hashref->{max_not_ordered_variants} if $hashref->{max_not_ordered_variants};
    $self->{max_not_ordered_variants_distance} = $hashref->{max_not_ordered_variants_distance} if $hashref->{max_not_ordered_variants_distance};
  }
  $self->{max_not_ordered_variants} = $Bio::EnsEMBL::VEP::Constants::MAX_NOT_ORDERED_VARIANTS if (!$self->{max_not_ordered_variants});
  $self->{max_not_ordered_variants_distance} = $Bio::EnsEMBL::VEP::Constants::MAX_NOT_ORDERED_VARIANTS_DISTANCE if (!$self->{max_not_ordered_variants_distance});

  return $self;
}


=head2 parser

  Arg 1      : (optional) Bio::EnsEMBL::VEP::Parser $parser
  Example    : $parser = $ib->parser();
  Description: Get/set the associated parser
  Returntype : Bio::EnsEMBL::VEP::Parser
  Exceptions : throws if given ref is not a Bio::EnsEMBL::VEP::Parser
  Caller     : next(), finish_annotation()
  Status     : Stable

=cut

sub parser {
  my $self = shift;

  if(@_) {
    my $parser = shift;
    assert_ref($parser, 'Bio::EnsEMBL::VEP::Parser');
    $self->{parser} = $parser;
  }

  return $self->{parser};
}


=head2 next

  Example    : $vfs = $ib->next();
  Description: Resets buffer, fills up to buffer_size from parser,
               returns reference to filled buffer. Returns emtpy
               arrayref at end of file.
  Returntype : Bio::EnsEMBL::VEP::Parser
  Exceptions : throws if given ref is not a Bio::EnsEMBL::VEP::Parser
  Caller     : next(), finish_annotation()
  Status     : Stable

=cut

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
  my $prev_start = 0;
  my $error_msg = "Exiting the program. The input file appears to be unsorted. Please sort by chromosome and by location and re-submit.\n";

  while(@$pre_buffer && @$buffer < $buffer_size) {
    my $vf = $pre_buffer->[0];

    # new chromosome
    if($prev_chr && $vf->{chr} ne $prev_chr) {
      $self->split_variants() if $self->{minimal};
      return $buffer;
    }

    # same chromosome
    else {
      push @$buffer, shift @$pre_buffer;
      $prev_chr = $vf->{chr};
    }
  }
  
  my $unsorted_formats = $self->config->{_params}->{unsorted_formats};
  
  if(my $parser = $self->parser) {
    while(@$buffer < $buffer_size && (my $vf = $parser->next)) {

      # skip long and unsupported types of SV; doing this here to avoid stopping looping
      next if $vf->{vep_skip};

      # exit the program if the maximum number of variants not ordered in the input file is reached
      if (!$self->param('no_check_variants_order') &&
          $self->{count_not_ordered_variants} &&
          $self->{count_not_ordered_variants} > $self->{max_not_ordered_variants}
      ) {
        die($error_msg);
      }

      # new chromosome
      if($prev_chr && $vf->{chr} ne $prev_chr) {

        # we can't push the VF back onto the parser, so add it to $pre_buffer
        # and it will get picked up on the following next() call
        push @$pre_buffer, $vf;
        
        $self->split_variants() if $self->{minimal};
        $prev_start = 0;
        return $buffer;
      }

      # same chromosome
      else {
        push @$buffer, $vf;
        $prev_chr = $vf->{chr};
        if (!$self->param('no_check_variants_order') && !$unsorted_formats->{$self->param('format')}) {
          # Use a default distance to check if the variant is still in the same region
          # even if it's not ordered with the previous variant location: we use the VF start and not the VCF start
          # (they can be different when the variant is a deletion for instance and/or when the alleles can be minimised)
          my $prev_start_with_range = $prev_start - $self->{max_not_ordered_variants_distance};
             $prev_start_with_range = 1 if ($prev_start_with_range < 0);
          if ($prev_start_with_range > $vf->{start}) {
            $self->{count_not_ordered_variants} ++;
          }
        }
        $prev_start = $vf->{start};
      }
    }
  }

  # exit the program if the maximum number of variants not ordered in the input file is reached (second point of exit)
  if (!$self->param('no_check_variants_order') &&
      $self->{count_not_ordered_variants} &&
      $self->{count_not_ordered_variants} > $self->{max_not_ordered_variants}
  ) {
    die($error_msg);
  }

  $self->split_variants() if $self->{minimal};

  return $buffer;
}


=head2 get_overlapping_vfs
  
  Arg 1      : int $start
  Arg 2      : int $end
  Example    : $vfs = $ib->get_overlapping_vfs($start, $end);
  Description: Gets all VariationFeatures overlapping the coordinate
               range given by $start, $end
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Caller     : AnnotationSource
  Status     : Stable

=cut

sub get_overlapping_vfs {
  my $self = shift;
  my $start = shift;
  my $end = shift;

  my @all_vfs;
  ($start, $end) = ($end, $start) if $start > $end;
  
  ## vfs obtained from trees are added to a separate array to allow for checking of overlaps
  ## from both shifted and unshifted positions
  
  if(my $tree = $self->interval_tree) {
    @all_vfs = @{$tree->fetch($start - 1, $end)};
  }
  else{
    $tree = $self->hash_tree unless($tree);  
    @all_vfs = values %{{
      map {$_->{_hash_tree_id} => $_}   # use _hash_tree_id to uniquify
      map {@{$tree->{$_} || []}} # tree might be empty
      (
        int($start / $HASH_TREE_SIZE) # start of range
        ..
        (int($end / $HASH_TREE_SIZE) + 1) # end of range
      )
    }};
  }      
  
  my @vfs;
  foreach my $vf (@all_vfs)  {
    if (overlap($vf->{start}, $vf->{end}, $start, $end)){
      push(@vfs, $vf);
    }
    if((defined($vf->{unshifted_end}) && (defined($vf->{unshifted_start})))) {
      if (overlap($vf->{unshifted_start}, $vf->{unshifted_end}, $start, $end)) {
        push(@vfs, $vf);
      }
    }
  }
  return [@vfs];
}


=head2 interval_tree
  
  Example    : $tree = $ib->interval_tree();
  Description: Gets interval tree populated with coordinates of
               VariationFeatures currently in the buffer.
  Returntype : Set::IntervalTree
  Exceptions : none
  Caller     : get_overlapping_vfs()
  Status     : Stable

=cut

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


=head2 hash_tree
  
  Example    : $tree = $ib->hash_tree();
  Description: Gets hash tree populated with coordinates of
               VariationFeatures currently in the buffer.

               This is a pure-perl interval tree of sorts;
               variants are added to "branches" (keys) of a hash,
               where the key is a bin obtained by
               
               bin = int(coord / $HASH_TREE_SIZE).
               
               Variants can span multiple bins and so can be added
               to more than one branch.
  Returntype : hashref
  Exceptions : none
  Caller     : get_overlapping_vfs()
  Status     : Stable

=cut

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


=head2 min_max
  
  Example    : my ($min, $max) = @{$ib->min_max()};
  Description: Gets the minimum and maximum coordinates spanned
               by variants in the buffer.
  Returntype : arrayref [$min, $max]
  Exceptions : none
  Caller     : AnnotationSource
  Status     : Stable

=cut

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


=head2 finish_annotation
  
  Example    : $ib->finish_annotation();
  Description: Finalises annotation on variants in the buffer, logs stats
               and ensures VariationFeatures have an attached slice.
  Returntype : none
  Exceptions : none
  Caller     : AnnotationSource
  Status     : Stable

=cut

sub finish_annotation {
  my $self = shift;
  $self->stats->log_lines_read($self->parser->line_number) if $self->parser;
  
  foreach my $vf(@{$self->buffer}) {
    $vf->{slice} ||= $self->get_slice($vf->{chr});
    $vf->_finish_annotation;
  }
}


=head2 reset_buffer
  
  Example    : $ib->reset_buffer();
  Description: Empties buffer and deletes interval/hash trees
  Returntype : none
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub reset_buffer {
  $_[0]->{temp} = {};
}


=head2 reset_pre_buffer
  
  Example    : $ib->reset_pre_buffer();
  Description: Empties pre-buffer
  Returntype : none
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub reset_pre_buffer {
  $_[0]->{pre_buffer} = [];
}


=head2 buffer
  
  Example    : $vfs = $ib->buffer();
  Description: Gets/sets arrayref of VariationFeatures that are the
               current buffer "fill"
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub buffer {
  my $self = shift;
  $self->{temp}->{buffer} = shift if @_;
  return $self->{temp}->{buffer} ||= [];
}


=head2 pre_buffer
  
  Example    : $vfs = $ib->pre_buffer();
  Description: Gets/sets arrayref of VariationFeatures that are the
               in the pre-buffer - used to fill buffer before reading
               from the parser
  Returntype : arrayref of Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub pre_buffer {
  return $_[0]->{pre_buffer} ||= [];
}


=head2 rejoin_required
  
  Arg 1      : (optional) bool $rejoin_required
  Example    : $ib->rejoin_required(1);
  Description: Gets/sets setting that informs downstream objects that
               variants in this buffer need re-joining after undergoing
               splitting in split_variants()
  Returntype : bool
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

sub rejoin_required {
  my $self = shift;
  $self->{temp}->{rejoin_required} = shift if @_;
  return $self->{temp}->{rejoin_required} ||= 0;
}


=head2 split_variants
  
  Example    : $ib->split_variants();
  Description: Splits VariationFeatures with multiple alternate alleles
               into one VariationFeature per alternate allele. Links are
               retained between split variants such that they can be
               rejoined later; this occurs in
               OutputFactory::rejoin_variants_in_InputBuffer().

               Splitting variants allows complex entries to be resolved
               to their minimum ALT/REF configuration before consequence
               calling e.g.:

               A/AC/G => -/C, A/G

               Note that coordinates are adjusted accordingly, meaning user
               must be careful interpreting output produced this way.
  Returntype : bool
  Exceptions : none
  Caller     : next()
  Status     : Stable

=cut

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

        ($ref, $alt, $start, $end, $this_changed) = @{trim_sequences($ref, $alt, $start, $end, 1)};
        $changed += $this_changed;

        # create a copy
        my $new_vf;
        %$new_vf = %{$original_vf};
        bless $new_vf, ref($original_vf);

        # give it a new allele string and coords
        $new_vf->allele_string($ref.'/'.$alt);
        $new_vf->{seq_region_start} = $new_vf->{start} = $start;
        $new_vf->{seq_region_end}   = $new_vf->{end}   = $end;
        $new_vf->{alt_allele}       = $alt;

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
