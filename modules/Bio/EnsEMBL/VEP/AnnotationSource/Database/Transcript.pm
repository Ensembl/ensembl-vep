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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript - database transcript annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript;

use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Database
  Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    core_type
    assembly
    gencode_basic
    all_refseq
    sift
    polyphen
    polyphen_analysis
    xref_refseq
    protein
    uniprot
    domains
  )]);

  $self->{cache_region_size} = 50000;
  $self->{source_type} = ($self->{core_type} || '') eq 'otherfeatures' ? 'refseq' : 'ensembl';

  return $self;
}

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;
  my @return;

  my $pfa = $self->get_adaptor('variation', 'PhenotypeFeature');

  my $cache_region_size = $self->{cache_region_size};

  foreach my $region(@{$regions}) {
    my ($c, $region_start) = @$region;

    my $slice = $self->get_slice($c);

    next unless $slice;

    # get a seq_region_Slice as for patch regions $slice won't cover the whole seq_region
    my $sr_slice = $slice->seq_region_Slice();

    my ($s, $e) = map {($_ - $slice->start) + 1} (
      ($region_start * $cache_region_size) + 1,
      ($region_start + 1) * $cache_region_size
    );

    # sanity check start and end
    $s = 1 if $s < 1;
    $e = $slice->length if $e > $slice->length;

    # get sub-slice
    my $sub_slice = $slice->sub_Slice($s, $e);

    next unless $sub_slice;

    # for some reason unless seq is called here the sequence becomes Ns later
    $sub_slice->seq;

    my @features;

    foreach my $gene(map {$_->transfer($sr_slice)} @{$sub_slice->get_all_Genes(undef, undef, 1)}) {
      my $gene_stable_id = $gene->stable_id;
      my $canonical_tr_id = $gene->{canonical_transcript_id};

      # any phenotypes?
      my $gene_has_phenotype = 0;

      if($pfa) {
        my $pfs = $pfa->fetch_all_by_Gene($gene);
        $gene_has_phenotype = $pfs && scalar @$pfs;
      }

      foreach my $tr(@{$gene->get_all_Transcripts}) {
        next unless $self->filter_transcript($tr);

        # in human and mouse otherfeatures DB, there may be duplicate genes
        # skip those from analysis refseq_human_import and refseq_mouse_import
        next if $self->{core_type} eq 'otherfeatures' && $self->{assembly} !~ /GRCh37/i && $tr->analysis && $tr->analysis->logic_name =~ /^refseq_[a-z]+_import$/;

        $tr->{_gene_stable_id} = $gene_stable_id;
        $tr->{_gene} = $gene;

        # indicate if canonical
        $tr->{is_canonical} = 1 if defined $canonical_tr_id and $tr->dbID eq $canonical_tr_id;

        # indicate phenotype
        $tr->{_gene_phenotype} = $gene_has_phenotype;

        $self->prefetch_transcript_data($tr);

        push @features, $tr;
      }
    }

    $cache->{$c}->{$region_start} = \@features;

    push @return, @features;
  }

  return \@return;
}

sub prefetch_transcript_data {
  my $self = shift;
  my $tr = shift;

  my $vep_cache = $tr->{_variation_effect_feature_cache} ||= {};

  $vep_cache->{introns} ||= $tr->get_all_Introns;
  $vep_cache->{sorted_exons} ||= [sort {$a->start <=> $b->start} @{$tr->get_all_Exons}];
  $vep_cache->{translateable_seq} ||= $tr->translateable_seq;
  $vep_cache->{mapper} ||= $tr->get_TranscriptMapper;

  # three prime UTR
  my $transferred = $tr->transfer($tr->feature_Slice());

  eval {
    $vep_cache->{three_prime_utr} = $transferred->three_prime_utr();
  };
  if($@) {
    warn "Problem getting 3' UTR:".$@;
  }

  # codon table
  unless ($vep_cache->{codon_table}) {
    # for mithocondrial dna we need to to use a different codon table
    my $attrib = $tr->slice->get_all_Attributes('codon_table')->[0];

    $vep_cache->{codon_table} = $attrib ? $attrib->value : 1;
  }

  # translation
  if(my $tl = $tr->translation) {
    $self->prefetch_translation_data($tr, $tl);
  }

  $self->prefetch_transcript_ids($tr);

  return $tr;
}

sub prefetch_translation_data {
  my $self = shift;
  my $tr = shift;
  my $tl = shift || $tr->translation;

  my $vep_cache = $tr->{_variation_effect_feature_cache} ||= {};

  # peptide
  unless ($vep_cache->{peptide}) {
    my $translation = $tr->translate;
    $vep_cache->{peptide} = $translation ? $translation->seq : undef;
  }

  # protein features
  if($self->{domains}) {
    my $pfs = $tr->translation ? $tr->translation->get_all_ProteinFeatures : [];

    # clean them to save cache space
    foreach my $pf(@$pfs) {

      # remove everything but the coord, analysis and ID fields
      foreach my $key(keys %$pf) {
        delete $pf->{$key} unless
          $key eq 'start' ||
          $key eq 'end' ||
          $key eq 'analysis' ||
          $key eq 'hseqname';
      }

      # remove everything from the analysis but the display label
      foreach my $key(keys %{$pf->{analysis}}) {
        delete $pf->{analysis}->{$key} unless $key eq '_display_label';
      }
    }

    $vep_cache->{protein_features} = $pfs;
  }

  # seq Edits
  $vep_cache->{seq_edits} = $tl->get_all_SeqEdits();

  # sift/polyphen
  my $pfpma = $self->get_adaptor('variation', 'ProteinFunctionPredictionMatrix');

  if($pfpma && defined($vep_cache->{peptide})) {
    foreach my $a('sift', 'polyphen_'.$self->{polyphen_analysis}) {
      next unless defined($self->{(split "_", $a)[0]});
      $vep_cache->{protein_function_predictions}->{$a} ||= $pfpma->fetch_by_analysis_translation_md5($a, md5_hex($vep_cache->{peptide}));
      delete $vep_cache->{protein_function_predictions}->{$a}->{adaptor};
    }
  }

  $self->prefetch_translation_ids($tr, $tl);

  return $tr;
}

sub prefetch_transcript_ids {
  my $self = shift;
  my $tr = shift;

  $tr->{_gene} ||= $tr->get_Gene();

  # gene symbol - get from gene cache if found already
  if(defined($tr->{_gene}->{_symbol})) {
    $tr->{_gene_symbol} = $tr->{_gene}->{_symbol};
    $tr->{_gene_symbol_source} = $tr->{_gene}->{_symbol_source};
    $tr->{_gene_hgnc_id} = $tr->{_gene}->{_hgnc_id}
  }
  else {
    $tr->{_gene_symbol} ||= undef;
    $tr->{_gene_symbol_source} ||= undef;

    if(my $xref = $tr->{_gene}->display_xref) {
      $tr->{_gene_symbol} = $xref->display_id;
      $tr->{_gene_symbol_source} = $xref->dbname;
      $tr->{_gene_hgnc_id} = $xref->primary_id if $xref->dbname eq 'HGNC';
    }

    else {
      my ($entry) = @{$tr->{_gene}->get_all_DBEntries('RefSeq_gene_name')};
      $tr->{_gene_symbol} = $entry->display_id if $entry;
    }

    # cache it on the gene object too
    $tr->{_gene}->{_symbol} = $tr->{_gene_symbol};
    $tr->{_gene}->{_symbol_source} = $tr->{_gene_symbol_source};
    $tr->{_gene}->{_hgnc_id} = $tr->{_gene_hgnc_id} if defined($tr->{_gene_hgnc_id});
  }

  # CCDS
  my @entries = grep {$_->database eq 'CCDS'} @{$tr->get_all_DBEntries};
  $tr->{_ccds} = $entries[0]->display_id if scalar @entries;
  $tr->{_ccds} ||= '-';

  # refseq
  if($self->{xref_refseq}) {
    @entries = grep {$_->database eq 'RefSeq_mRNA'} @{$tr->get_all_DBEntries};
    if(scalar @entries) {
      $tr->{_refseq} = join ",", map {$_->display_id} @entries;
    }
    else {
      $tr->{_refseq} = '-';
    }
  }
}

sub prefetch_translation_ids {
  my $self = shift;
  my $tr = shift;
  my $tl = shift || $tr->translation;

  # uniprot
  if($self->{uniprot}) {
    
    $tr->{_swissprot} = '-';
    my @entries = grep {$_->database eq 'Uniprot/SWISSPROT'} @{$tl->get_all_DBEntries};
    if(scalar @entries) {
      $tr->{_swissprot} = join ",", map {$_->display_id} @entries;
    }

    $tr->{_trembl} = '-';
    @entries = grep {$_->database eq 'Uniprot/SPTREMBL'} @{$tl->get_all_DBEntries};
    if(scalar @entries) {
      $tr->{_trembl} = join ",", map {$_->display_id} @entries;
    }


    $tr->{_uniparc} = '-';
    @entries = grep {$_->database eq 'UniParc'} @{$tl->get_all_DBEntries};
    if(scalar @entries) {
      $tr->{_uniparc} = join ",", map {$_->display_id} @entries;
    }
  }

  # Ensembl protein ID
  if($self->{protein}) {
    $tr->{_protein} = $tl->stable_id;
  }

  return $tr;
}

1;