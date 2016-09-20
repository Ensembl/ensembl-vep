=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript - base class for transcript annotation sources

DO NOT USE DIRECTLY

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::VEP::AnnotationSource);

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my $up_size   = $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE;
  my $down_size = $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE;
  my $tva = $self->get_adaptor('variation', 'TranscriptVariation');

  foreach my $tr(@{$self->get_all_features_by_InputBuffer($buffer)}) {
    my $tr_strand = $tr->{strand} || $tr->strand;
    my $fs = $tr->{start} - ($tr_strand == 1 ? $up_size : $down_size);
    my $fe = $tr->{end} + ($tr_strand == 1 ? $down_size : $up_size);
    my $slice;

    my $vfs = $buffer->get_overlapping_vfs($fs, $fe);
    next unless @$vfs;
    $tr = $self->lazy_load_transcript($tr);
    next unless $tr;

    foreach my $vf(@$vfs) {
      $vf->{slice} ||= $slice ||=  $tr->{slice};

      if(ref($vf) eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature') {
        my $svo = Bio::EnsEMBL::Variation::TranscriptStructuralVariation->new(
          -transcript                   => $tr,
          -structural_variation_feature => $vf,
          -no_transfer                  => 1
        );

        $vf->add_TranscriptStructuralVariation($svo);
      }

      else {
        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
          -transcript        => $tr,
          -variation_feature => $vf,
          -adaptor           => $tva,
          -no_ref_check      => 1,
          -no_transfer       => 1
        );

        $vf->add_TranscriptVariation($tv);
      }
    }
  }
}

sub up_down_size {
  return $_[0]->{up_down_size} ||= (
    sort {$a <=> $b} (
      $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE,
      $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE
    )
  )[-1];
}

sub filter_transcript {
  my $self = shift;
  my $t = shift;

  # there are some transcripts in the otherfeatures DB with no stable ID!!!
  return 0 unless $t->stable_id;

  # using gencode basic?
  if($self->{gencode_basic} && !(grep {$_->{code} eq 'gencode_basic'} @{$t->get_all_Attributes})) {
    return 0;
  }

  # using all_refseq?
  if(
    !$self->{all_refseq} &&
    (
      # we only want RefSeq transcripts e.g. NM_12930
      (
        $self->{source_type} eq 'refseq' &&
        ($t->stable_id || '') !~ /^[A-Z]{2}\_\d+/
      ) ||

      # and the same from the merged cache
      (
        $self->{source_type} eq 'merged' &&
        ($t->{_source_cache} || '') eq 'RefSeq' &&
        ($t->stable_id || '') !~ /^[A-Z]{2}\_\d+/
      )
    )
  ) {
    return 0;
  }

  return 1;
}
  
sub merge_features {
  my $self = shift;
  my $features = shift;

  my %hgnc_ids = ();
  my %refseq_stuff = ();
  my %seen_trs;
  my @return;

  my $source_type_is_refseq = $self->{source_type} && ($self->{source_type} eq 'refseq' || $self->{source_type} eq 'merged') ? 1 : 0;

  while(my $tr = shift @$features) {

    # there are some transcripts in the otherfeatures DB with no stable ID!
    next unless $tr->stable_id;

    # track already added transcripts by dbID
    my $dbID = $tr->dbID;
    if($seen_trs{$dbID}) {
      # $count_duplicates++;
      next;
    }

    ## hack to copy HGNC IDs
    $hgnc_ids{$tr->{_gene_symbol}} = $tr->{_gene_hgnc_id} if defined($tr->{_gene_hgnc_id});

    ## hack to copy RefSeq gene stuff
    if($source_type_is_refseq) {
      $refseq_stuff{$tr->{_gene}->stable_id}->{$_} ||= $tr->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);
    }

    $seen_trs{$dbID} = 1;

    push @return, $tr;
  }

  ## hack to copy HGNC IDs and RefSeq stuff
  foreach my $tr(@return) {
    $tr->{_gene_hgnc_id} = $hgnc_ids{$tr->{_gene_symbol}} if defined($tr->{_gene_symbol}) && defined($hgnc_ids{$tr->{_gene_symbol}});

    if($source_type_is_refseq) {
      $tr->{$_} ||= $refseq_stuff{$tr->{_gene}->stable_id}->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);
    }
  }

  return \@return;
}

sub lazy_load_transcript {
  my ($self, $tr) = @_;
  if($tr->{_vep_lazy_loaded}) {
    return $tr;
  }
  else {
    $tr->{_vep_lazy_loaded} = 1;
    return $tr;
  }
}

1;