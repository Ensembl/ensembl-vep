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

  if(ref($as) =~ /Database/) {
    my $ta = $self->get_adaptor('core', 'transcript');
    $self->insert($_->seq_region_name, $_->seq_region_start, $_->seq_region_end) for @{$ta->fetch_all_by_biotype('protein_coding')};
  }
  elsif(ref($as) =~ /Cache/) {
    my $as_dir = $as->dir;

    if(-e $as_dir.'/transcript_coords.txt') {
      open TR, $as_dir.'/transcript_coords.txt';
      while(<TR>) {
        chomp;
        $self->insert(split);
      }
      close TR;
    }
    else {
      my $cache_region_size = $as->{cache_region_size};

      open TR, ">".$as_dir.'/transcript_coords.txt';

      opendir DIR, $as_dir;
      foreach my $c(grep {!/^\./ && -d $as_dir.'/'.$_} readdir DIR) {
        
        opendir CHR, $as_dir.'/'.$c;
        foreach my $file(grep {/\d+\-\d+\.gz/} readdir CHR) {
          my ($s) = split(/\D/, $file);

          foreach my $t(
            grep {$_->biotype eq 'protein_coding'}
            @{$as->get_features_by_regions_uncached([[$c, ($s - 1) / $cache_region_size]])}
          ) {
            my ($s, $e) = ($t->seq_region_start, $t->seq_region_end);
            $self->insert($c, $s, $e);
            print TR "$c\t$s\t$e\n";
          }
        }
        closedir CHR;
      }
      closedir DIR;
      close TR;
    }
  }
  elsif(ref($as) =~ /File/) {
    my $parser = $as->parser;

    foreach my $chr(@{$parser->{tabix_file}->seqnames}) {
      $parser->seek($chr, 1, 1e10);
      while($parser->next) {
        $self->insert($chr, $parser->get_start, $parser->get_end);
      }
    }

    delete $as->{parser};
  }  
  else {
    throw("ERROR: Don't know how to process ".ref($as));
  }

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
  return $self->get_chr_tree($c)->fetch($s - 1, $e);
}

1;