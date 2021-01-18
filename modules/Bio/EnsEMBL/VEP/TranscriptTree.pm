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

# EnsEMBL module for Bio::EnsEMBL::VEP::TranscriptTree
#
#

=head1 NAME

Bio::EnsEMBL::VEP::TranscriptTree - class containing IntervalTree of transcript locations

=head1 SYNOPSIS

my $tree = Bio::EnsEMBL::VEP::TranscriptTree->new({
  config            => $config,
  annotation_source => $as
});

my $overlaps = $tree->fetch($chr, $start, $end);

=head1 DESCRIPTION

The TranscriptTree class uses a series of Set::IntervalTree objects
to retrieve references to overlapping features.

It requires an attached AnnotationSource that implements a
populate_tree() method that will populate the tree with features.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::TranscriptTree;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Set::IntervalTree;


=head2 new

  Arg 1      : hashref $args
               {
                 config            => Bio::EnsEMBL::VEP::Config,
                 annotation_source => Bio::EnsEMBL::VEP::AnnotationSource,
               }
  Example    : $tree = Bio::EnsEMBL::VEP::TranscriptTree->new({
                 config            => $config,
                 annotation_source => $as
               });
  Description: Create a new Bio::EnsEMBL::VEP::TranscriptTree object. Requires
               annotation_source arg to be set to a reference to a compatible
               AnnotationSource i.e. one that implements a populate_tree()
               method.
  Returntype : Bio::EnsEMBL::VEP::TranscriptTree
  Exceptions : throws if no AnnotationSource given or AnnotationSource has no populate_tree() method
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  my $as = $hashref->{annotation_source};
  assert_ref($as, 'Bio::EnsEMBL::VEP::AnnotationType::Transcript');

  throw("ERROR: Unable to populate tree from annotation source type ".ref($as)) unless $as->can('populate_tree');
  $as->populate_tree($self);

  $self->valid_chromosomes($as->valid_chromosomes);

  return $self;
}


=head2 get_chr_tree

  Example    : $chr_tree = $tree->get_chr_tree($chr);
  Description: Get the Set::IntervalTree object for this chromosome.
  Returntype : Set::IntervalTree
  Exceptions : none
  Caller     : insert(), fetch()
  Status     : Stable

=cut

sub get_chr_tree {
  my $self = shift;
  my $chr = shift;
  return $self->{trees}->{$chr} ||= Set::IntervalTree->new();
}


=head2 insert

  Arg 1      : string $chromosome
  Arg 2      : int $start
  Arg 3      : int $end
  Arg 4      : scalar $obj
  Example    : [$min, $max] = $tree->insert($chr, $start, $end);
               $tree->insert($chr, $start, $end, $obj);
  Description: Insert given chromosome range into tree. If no object $obj
               is provided, a fetch is first performed, and any overlapping
               ranges are merged with this one, removed, and a single
               expanded range re-inserted into the tree. The expanded range
               is then returned as an arrayref [$min, $max].

               If an object $obj (can be a reference or any scalar) then this
               is inserted instead, with the provided [$start, $end] returned
  Returntype : arrayref [$min, $max]
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub insert {  
  my ($self, $c, $s, $e, $obj) = @_;

  my $tree = $self->get_chr_tree($c);

  # obj given, just insert it and return
  if($obj) {
    $obj->{s} ||= $s;
    $obj->{e} ||= $e;

    $tree->insert($obj, $s - 1, $e);
    return [$s, $e];
  }

  # otherwise do a merge if any overlap and insert an arrayref representing the range
  my ($min, $max) = ($s, $e);

  my $fetched = $tree->fetch($s - 1, $e);
  if(@$fetched) {
    foreach my $f(@$fetched) {
      $min = $f->[0] if $f->[0] < $min;
      $max = $f->[1] if $f->[1] > $max;
    }

    $tree->remove($s - 1, $e);
  }

  $tree->insert([$min, $max], $min - 1, $max);

  return [$min, $max];
}


=head2 fetch

  Arg 1      : string $chromosome
  Arg 2      : int $start
  Arg 3      : int $end
  Example    : my $overlapping = $tree->fetch($chr, $start, $end)
  Description: Fetch features overlapping the given coordinate range
  Returntype : arrayref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch {
  my ($self, $c, $s, $e) = @_;
  return $self->get_chr_tree($self->get_source_chr_name($c))->fetch($s - 1, $e);
}


=head2 nearest

  Arg 1      : string $chromosome
  Arg 2      : int $start
  Arg 3      : int $end
  Example    : my $nearest = $tree->nearest($chr, $start, $end)
  Description: Fetch nearest features to given coordinate range
  Returntype : whatever reference type was inserted by insert()
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub nearest {
  my ($self, $c, $s, $e) = @_;

  my $tree = $self->get_chr_tree($self->get_source_chr_name($c));

  # invert if $s > $e
  ($e, $s) = ($s, $e) if $s > $e;

  # search up and down - note "up" in Set::IntervalTree means a higher numerical value coordinate
  # confusing if you're used to talking in terms of up/downstream on a genome
  # so we search "up" from $s and "down" from $e as this will capture overlaps too
  my @search = grep {$_} (
    $tree->fetch_nearest_down($e),
    $tree->fetch_nearest_up($s),
  );

  # bail out here if <=1 result (empty chr or 1 result)
  return \@search unless scalar @search > 1;

  # get dists as hash indexed on array pos in @search
  my %dists =
    map {$_ => $self->_get_dist(
      @{$self->_get_obj_start_end($search[$_])},
      $s, $e
    )}
    0..$#search;

  # find minimum distance
  my $min_dist = (sort {$a <=> $b} values %dists)[0];

  # return all objects that have that dist
  return [map {$search[$_]} grep {$dists{$_} == $min_dist} 0..$#search];
}


=head2 _get_obj_start_end

  Arg 1      : hashref $obj or arrayref [$start, $end]
  Example    : [$s, $e] = $tree->_get_obj_start_end($obj)
  Description: Gets start and end given either an object hashref
               or arrayref [$s, $e]
  Returntype : arrayref [$start, $end]
  Exceptions : none
  Caller     : nearest()
  Status     : Stable

=cut

sub _get_obj_start_end {
  my ($self, $obj) = @_;

  my ($s, $e);
  
  if(ref($obj) eq 'ARRAY') {
    ($s, $e) = @$obj;
  }
  else {
    $s = $obj->{s} || $obj->{start} or throw("ERROR: No start field \"s\" defined on object\n");
    $e = $obj->{e} || $obj->{end} || $s;
  }

  return [$s, $e];
}


=head2 _get_dist

  Arg 1      : int $start1
  Arg 2      : int $end1
  Arg 3      : int $start2
  Arg 4      : int $end2
  Example    : $dist = $tree->_get_dist(1, 5, 8, 12)
  Description: Finds shortest distance between any two ends of two given
               coordinate pairs.
  Returntype : int
  Exceptions : none
  Caller     : nearest()
  Status     : Stable

=cut

sub _get_dist {
  my ($self, $o_s, $o_e, $s, $e) = @_;

  # calculate all possible distances
  # sort to find lowest
  return (
    sort {$a <=> $b}
    (
      abs($s - $o_s),
      abs($s - $o_e),
      abs($e - $o_s),
      abs($e - $o_e)
    )
  )[0];
}


=head2 valid_chromosomes
  
  Arg 1      : (optional) arrayref $valid_chromosomes
  Example    : $valids = $tree->valid_chromosomes();
  Description: Getter/setter for the list of valid chromosomes as found
               in the configured AnnotationSource.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : new(), BaseVEP
  Status     : Stable

=cut

sub valid_chromosomes {
  my $self = shift;
  $self->{valid_chromosomes} = shift if @_;
  return $self->{valid_chromosomes};
}

1;