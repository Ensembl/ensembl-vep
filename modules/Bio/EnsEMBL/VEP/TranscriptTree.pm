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

# EnsEMBL module for Bio::EnsEMBL::VEP::TranscriptTree
#
#

=head1 NAME

Bio::EnsEMBL::VEP::TranscriptTree - class containing IntervalTree of transcript locations

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::TranscriptTree;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Set::IntervalTree;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  my $as = $hashref->{annotation_source};
  assert_ref($as, 'Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript');

  throw("ERROR: Unable to populate tree from annotation source type ".ref($as)) unless $as->can('populate_tree');
  $as->populate_tree($self);

  $self->valid_chromosomes($as->valid_chromosomes);

  return $self;
}

sub get_chr_tree {
  my $self = shift;
  my $chr = shift;
  return $self->{trees}->{$chr} ||= Set::IntervalTree->new();
}

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

sub fetch {
  my ($self, $c, $s, $e) = @_;
  return $self->get_chr_tree($self->get_source_chr_name($c))->fetch($s - 1, $e);
}

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

sub valid_chromosomes {
  my $self = shift;
  $self->{valid_chromosomes} = shift if @_;
  return $self->{valid_chromosomes};
}

1;