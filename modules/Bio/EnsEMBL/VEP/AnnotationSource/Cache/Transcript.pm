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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript - local disk transcript annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Variation::TranscriptVariation;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCache);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(gencode_basic all_refseq)]); 

  return $self;
}

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my $up_size   = $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE;
  my $down_size = $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE;
  my $tva = $self->get_adaptor('variation', 'TranscriptVariation');

  foreach my $tr(@{$self->get_all_features_by_InputBuffer($buffer, $self->{cache_region_size})}) {
    my $fs = $tr->{start} - ($tr->strand == 1 ? $up_size : $down_size);
    my $fe = $tr->{end} + ($tr->strand == 1 ? $down_size : $up_size);

    foreach my $vf(grep { overlap($fs, $fe, $_->{start}, $_->{end}) } @{$buffer->buffer}) {
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

sub get_dump_file_name {
  my $self = shift;
  my $chr  = shift;
  my $region = shift;

  throw("No chromosome given") unless $chr;
  throw("No region given") unless $region;

  # allow to pass region (start-end) or $start, $end
  $region .= '-'.shift if @_;

  return sprintf(
    "%s/%s/%s\.%s",
    $self->dir,
    $chr,
    $region,
    $self->file_suffix
  );
}

sub deserialized_obj_to_features {
  my $self = shift;
  my $obj = shift;

  my $tra = $self->get_adaptor('core', 'Translation');
  my $sa  = $self->get_adaptor('core', 'Slice');

  my @features;

  foreach my $chr(keys %$obj) {
    foreach my $t(@{$obj->{$chr}}) {

      # we may need to filter some out based on config options
      # do this ASAP to avoid processing more features than we need to
      next unless $self->filter_transcript($t);

      if(defined($t->{translation})) {
        $t->{translation}->{adaptor} = $tra;
        $t->{translation}->{transcript} = $t;
        weaken($t->{translation}->{transcript});
      }

      $t->{slice}->{adaptor} = $sa;

      $_->{slice} ||= $t->{slice} for @{$t->{_trans_exon_array}};

      push @features, $t;
    }
  }

  return \@features;
}

sub filter_transcript {
  my $self = shift;
  my $t = shift;

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
    if($self->{source_type} eq 'refseq' || $self->{source_type} eq 'merged') {
      $refseq_stuff{$tr->{_gene}->stable_id}->{$_} ||= $tr->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);
    }

    $seen_trs{$dbID} = 1;

    push @return, $tr;
  }

  ## hack to copy HGNC IDs and RefSeq stuff
  foreach my $tr(@return) {
    $tr->{_gene_hgnc_id} = $hgnc_ids{$tr->{_gene_symbol}} if defined($tr->{_gene_symbol}) && defined($hgnc_ids{$tr->{_gene_symbol}});

    if($self->{source_type} eq 'refseq' || $self->{source_type} eq 'merged') {
      $tr->{$_} ||= $refseq_stuff{$tr->{_gene}->stable_id}->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);
    }
  }

  return \@return;
}

1;