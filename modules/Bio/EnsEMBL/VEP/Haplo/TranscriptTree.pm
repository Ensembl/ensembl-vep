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

# EnsEMBL module for Bio::EnsEMBL::VEP::Haplo::TranscriptTree
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Haplo::TranscriptTree - class containing IntervalTree of transcript locations

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Haplo::TranscriptTree;

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

  $self->valid_chromosomes({map {$_ => 1} @{$as->get_valid_chromosomes}});

  return $self;
}

sub get_chr_tree {
  my $self = shift;
  my $chr = shift;
  return $self->{trees}->{$chr} ||= Set::IntervalTree->new();
}

sub insert {  
  my ($self, $c, $s, $e) = @_;

  my $tree = $self->get_chr_tree($c);

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

sub valid_chromosomes {
  my $self = shift;
  $self->{valid_chromosomes} = shift if @_;
  return $self->{valid_chromosomes};
}

sub get_source_chr_name {
  my $self = shift;
  my $chr = shift;

  my $chr_name_map = $self->{_chr_name_map} ||= {};

  if(!exists($chr_name_map->{$chr})) {
    my $mapped_name = $chr;

    my $valid = $self->valid_chromosomes;

    unless($valid->{$chr}) {

      # try synonyms first
      my $synonyms = $self->chromosome_synonyms;

      foreach my $syn(keys %{$synonyms->{$chr} || {}}) {
        if($valid->{$syn}) {
          $mapped_name = $syn;
          last;
        }
      }

      # still haven't got it
      if($mapped_name eq $chr) {

        # try adding/removing "chr"
        if($chr =~ /^chr/i) {
          my $tmp = $chr;
          $tmp =~ s/^chr//i;

          $mapped_name = $tmp if $valid->{$tmp};
        }
        elsif($valid->{'chr'.$chr}) {
          $mapped_name = 'chr'.$chr;
        }
      }
    }

    $chr_name_map->{$chr} = $mapped_name;
  }

  return $chr_name_map->{$chr};
}

1;