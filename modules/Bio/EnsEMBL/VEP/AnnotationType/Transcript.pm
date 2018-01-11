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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationType::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationType::Transcript - base class for transcript annotation sources

=head1 SYNOPSIS

Should not be invoked directly.

=head1 DESCRIPTION

Helper class for all transcript-based AnnotationSource classes. Contains the bulk of the
code that carries out the annotation.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationType::Transcript;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);


our ($CAN_USE_HTS, $CAN_USE_INTERVAL_TREE);

BEGIN {
  if (eval q{ require Bio::DB::HTS; 1 }) {
    $CAN_USE_HTS = 1;
  }
  if (eval q{ require Bio::EnsEMBL::VEP::TranscriptTree; 1 }) {
    $CAN_USE_INTERVAL_TREE = 1;
  }
}

our $DEBUG = 0;


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $as->annotate_InputBuffer($ib);
  Description: Gets overlapping transcripts for the variants in
               the input buffer, and creates TranscriptVariation objects for
               each overlapping pair of variant/transcript. These are then added
               to the relevant VariationFeature.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my $up_size   = $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE;
  my $down_size = $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE;
  my $tva = $self->get_adaptor('variation', 'TranscriptVariation');
  my $use_feature_ref = $self->{use_transcript_ref};

  foreach my $tr(@{$self->get_all_features_by_InputBuffer($buffer)}) {
    my $tr_strand = $tr->{strand} || $tr->strand;
    my $fs = $tr->{start} - ($tr_strand == 1 ? $up_size : $down_size);
    my $fe = $tr->{end} + ($tr_strand == 1 ? $down_size : $up_size);
    my $slice;

    # get overlapping VFs
    my $vfs = $buffer->get_overlapping_vfs($fs, $fe);
    next unless @$vfs;

    # lazy load transcript
    $tr = $self->lazy_load_transcript($tr);
    next unless $tr;

    # apply filters
    next if $self->filter_set && !$self->filter_set->evaluate($tr);

    # apply edits
    $self->apply_edits($tr);

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
          -no_transfer       => 1,
          -use_feature_ref   => $use_feature_ref,
        );

        $vf->add_TranscriptVariation($tv);
      }
    }
  }

  if(my $nearest_type = $self->{nearest}) {
    foreach my $vf(@{$buffer->buffer}) {
      $vf->{nearest} = $self->get_nearest($vf, $nearest_type);
    }
  }
}


=head2 up_down_size

  Example    : $size = $as->up_down_size();
  Description: Gets range in bp that should be added to boundaries
               when fetching features.
  Returntype : int
  Exceptions : none
  Caller     : get_all_regions_by_InputBuffer()
  Status     : Stable

=cut

sub up_down_size {
  return $_[0]->{up_down_size} ||= (
    sort {$a <=> $b} (
      $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE,
      $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE
    )
  )[-1];
}


=head2 filter_transcript

  Arg 1      : Bio::EnsEMBL::Transcript $t
  Example    : $pass = $as->filter_transcript($t);
  Description: Pass/fail a transcript based on user configuration.
               Transcripts may be excluded if they are not in the GENCODE
               basic set, or if they are not RefSeq transcripts from the
               refseq cache.
  Returntype : bool
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

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
      # we only want RefSeq transcripts e.g. NM_12930,
      # or 4540 for MT transcripts
      (
        $self->{source_type} eq 'refseq' &&
        ($t->stable_id || '') !~ /^[A-Z]{2}\_\d+|^\d{4}$/
      ) ||

      # and the same from the merged cache
      (
        $self->{source_type} eq 'merged' &&
        ($t->{_source_cache} || '') eq 'RefSeq' &&
        ($t->stable_id || '') !~ /^[A-Z]{2}\_\d+|^\d{4}$/
      )
    )
  ) {
    return 0;
  }

  return 1;
}


=head2 merge_features

  Arg 1      : arrayref of Bio::EnsEMBL::Transcript $transcripts
  Example    : $merged = $as->merge_features($features);
  Description: "Uniquifies" a list of transcripts using dbID
  Returntype : arrayref of Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut
  
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
  my %by_stable_id;
  foreach my $tr(@return) {
    $tr->{_gene_hgnc_id} = $hgnc_ids{$tr->{_gene_symbol}} if defined($tr->{_gene_symbol}) && defined($hgnc_ids{$tr->{_gene_symbol}});

    if($source_type_is_refseq) {
      $tr->{$_} ||= $refseq_stuff{$tr->{_gene}->stable_id}->{$_} for qw(_gene_symbol _gene_symbol_source _gene_hgnc_id);
      push @{$by_stable_id{$tr->{stable_id}}}, $tr;
    }
  }

  ## now remove duplicates...
  if($source_type_is_refseq) {
    my @new;
    my %done_stable_id = ();

    foreach my $tr(@return) {
      my $stable_id = $tr->{stable_id};
      next if $done_stable_id{$stable_id};
      $done_stable_id{$stable_id} = 1;

      my @all_with_stable_id = @{$by_stable_id{$stable_id}};

      if(scalar @all_with_stable_id == 1) {
        push @new, @all_with_stable_id;
      }
      else {
        # try and find one with source 'ensembl'
        my ($ensembl_tr) = grep {$_->{source} eq 'ensembl'} @all_with_stable_id;

        if($ensembl_tr) {
          push @new, $ensembl_tr;
        }
        else {
          # just take the one with the lowest dbID
          push @new, (sort {$a->dbID <=> $b->dbID} @all_with_stable_id)[0];
        }
      }
    }

    @return = @new;
  }

  return \@return;
}


=head2 lazy_load_transcript

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Example    : $as->lazy_load_transcript($tr);
  Description: Stub method to lazy-load transcript data
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

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



## BAM EDIT METHODS
###################


=head2 bam

  Arg 1      : (optional) string $bam_file
  Example    : $bam_obj = $as->bam($file);
  Description: Get Bio::DB::HTS object for a BAM file, used for
               editing transcripts in apply_edits()
  Returntype : Bio::DB::HTS
  Exceptions : throws if Bio::DB::HTS not installed
  Caller     : apply_edits()
  Status     : Stable

=cut

sub bam {
  my $self = shift;

  if(@_ || (!$self->{_bam} && $self->{bam})) {
    throw("ERROR: Cannot add BAM file without Bio::DB::HTS installed\n") unless $CAN_USE_HTS;
    $self->{_bam} = Bio::DB::HTS->new(-bam => @_ ? shift : $self->{bam});
  }

  $self->{_bam} = undef unless exists $self->{_bam};

  return $self->{_bam};
}


=head2 apply_edits

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Arg 2      : (optional) Bio::DB::HTS $bam
  Example    : $as->apply_edits($tr);
  Description: Uses a BAM file containing aligned sequence for transcripts
               to "correct" the sequence retrieved for a transcript from the
               underlying genome. Edits are applied as SeqEdit attributes,
               which then get picked up when calling e.g. $tr->spliced_seq
  Returntype : Bio::EnsEMBL::Transcript unless fails
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub apply_edits {
  my ($self, $tr, $bam) = @_;

  return $tr if exists($tr->{_bam_edit_status});

  $bam ||= $self->bam;
  return unless $bam;

  my $stable_id = $tr->stable_id;

  # COMMENTED OUT FOR NOW IN CASE WE CAN'T TRUST THESE ATTRIBS
  # don't need to edit if exact match already
  if(grep {$_->code eq 'rseq_mrna_match'} @{$tr->get_all_Attributes}) {
    print STDERR "OK $stable_id MATCH ATTRIBUTE PRESENT\n" if $DEBUG;
    return $tr;
  }

  # get the alignment representing this transcript aligned to the genome
  my ($al) =
    grep {$_->query->name eq $stable_id}
    ($bam->get_features_by_location(
      -seq_id => $self->get_source_chr_name($tr->seq_region_name, 'bam', [$bam->seq_ids]),
      -start  => $tr->start,
      -end    => $tr->end,
      -types  => 'match',
    ));

  unless($al) {
    print STDERR "FAILED $stable_id NO ALIGNMENT\n" if $DEBUG;
    return;
  }

  # get cigar string and check for indels and mismatches
  my $cigar_array    = $al->cigar_array;
  my $q_seq_obj      = $al->query->seq;
  my $mapping_strand = $al->strand;

  # we need to look at these op types
  my %edit_ops = map {$_ => 1} qw(X D I);

  ## COMMENTED OUT FOR NOW, TRUST SEQ COMPARISON INSTEAD
  ## check for edit ops
  # unless(grep {$edit_ops{$_->[0]}} @$cigar_array) {
  #   print STDERR "OK $stable_id NO EDITS REQUIRED\n" if $DEBUG;
  #   return;
  # }

  # do a manual check on the sequences
  my $pre_edit_seq = $tr->spliced_seq;
  my $bam_seq = ($mapping_strand > 0 ? $q_seq_obj->seq : $q_seq_obj->revcom->seq);

  if($pre_edit_seq eq $bam_seq) {
    print STDERR "OK $stable_id SEQUENCES MATCH ANYWAY\n" if $DEBUG;
    return $tr;
  }

  # op_consumes is a hash saying whether each CIGAR op type "consumes" query [0] and/or ref [1] seq
  # e.g. M (match) => [1, 1], I (insertion) => [1, 0], D (deletion) => [0, 1]
  # modified from Bio::Cigar
  my %op_consumes = (
    # op => [query, reference]
    'M' => [1, 1],
    'I' => [1, 0],
    'D' => [0, 1],
    'N' => [0, 0], # special case intron - normally this should be [0, 1] but the position we want is cDNA-relative
    'S' => [1, 0],
    'H' => [0, 0],
    'P' => [0, 0],
    '=' => [1, 1],
    'X' => [1, 1],
  );

  my @edits;
  my $q_seq    = $q_seq_obj->seq;
  my $bam_file = (split('/', $self->{bam}))[-1];
  my %seen_ops = ();

  # work out positions to start from
  my $current_q_pos = 1;
  my $current_t_pos;

  # forward strand
  if($mapping_strand > 0) {
    $current_t_pos = ($al->start - $tr->seq_region_start) + 1;
  }
  # reverse strand
  else {
    $current_t_pos = ($al->start - $tr->seq_region_start) + length($tr->spliced_seq);
  }

  foreach my $c(@$cigar_array) {
    my ($op, $l) = @$c;

    throw("ERROR: Unrecognised operation $op\n") unless $op_consumes{$op};

    if($edit_ops{$op}) {
      $seen_ops{$op}++;

      my $q_s = $current_q_pos;

      my ($t_s, $t_e);

      if($op eq 'I') {
        if($mapping_strand > 0) {
          $t_s = $current_t_pos;
          $t_e = $current_t_pos - 1;
        }
        else {
          $t_s = $current_t_pos + 1;
          $t_e = $current_t_pos;
        }
      }
      else {
        if($mapping_strand > 0) {
          $t_s = $current_t_pos;
          $t_e = ($current_t_pos + $l) - 1;
        }
        else {
          $t_s = ($current_t_pos - $l) + 1;
          $t_e = $current_t_pos;
        }
      }

      # get query seq, account for dels by multiplying length by the op_consumes value for query seq
      my $alt_seq = substr($q_seq, $current_q_pos - 1, $op_consumes{$op}->[0] * $l);

      # reverse complement if necessary
      reverse_comp(\$alt_seq) if $mapping_strand < 0;

      push @edits, Bio::EnsEMBL::Attribute->new(
        -VALUE       => "$t_s $t_e $alt_seq",
        -CODE        => '_rna_edit',
        -NAME        => 'RNA Edit',
        -DESCRIPTION => "Edit from $bam_file, op=$op len=$l mapped to $t_s-$t_e"
      );
    }

    # increment positions
    $current_q_pos += $l if $op_consumes{$op}->[0];
    $current_t_pos += ($l * $mapping_strand) if $op_consumes{$op}->[1];
  }

  # add the edits to the transcript as attributes
  $tr->add_Attributes(@edits);

  # tell the transcript to apply them when we call relevant methods
  my $edits_enabled_bak = $tr->edits_enabled(); 
  $tr->edits_enabled(1);

  # now test whether the sequence the transcript object will create with the SeqEdits matches the input (RefSeq) seq
  my $new_tr_spliced_seq = '';
  eval {$new_tr_spliced_seq = $tr->spliced_seq};

  my $cmp = !$@ && $new_tr_spliced_seq eq $bam_seq;
  # if($DEBUG) {
  #   my $status = $cmp ? 'OK' : 'FAILED';
  #   my $error = $@ ? "\t$@" : "";
  #   $error =~ s/\s+$//g;

  #   print STDERR
  #     "$status $stable_id STRAND $mapping_strand OPS ".
  #     join(", ", map {$_.":".$seen_ops{$_}} sort keys %seen_ops).
  #     " EDITS ".
  #     (join(", ", map {$_->value} grep {$_->code eq '_rna_edit'} @{$tr->get_all_Attributes}) || 'NONE').
  #     "$error\n";

  #   if($status eq 'FAILED') {
  #     open OUT, ">$stable_id.fa";
  #     print OUT "\>BAM\n$bam_seq\n\>PRE\n$pre_edit_seq\n\>TR\n".$new_tr_spliced_seq;
  #     close OUT;
  #   }
  # }
  
  if($cmp) {
    # flag successful
    $tr->{_bam_edit_status} = 'ok';

    # delete stuff pre-cached for VEP as they need to be regenerated with edited seq
    if(my $vep_cache = $tr->{_variation_effect_feature_cache}) {
      delete $vep_cache->{$_} for qw(
        mapper
        five_prime_utr
        three_prime_utr
        translateable_seq
        peptide
        protein_function_predictions
        protein_features
      );
    }

    # now add new spliced_seq field as this will be needed by use_transcript_ref
    $tr->{_variation_effect_feature_cache}->{spliced_seq} = $new_tr_spliced_seq;
  }
  else {
    # flag as failed
    $tr->{_bam_edit_status} = 'failed';

    # revert the changes we've made and flag this transcript
    $tr->{attributes} = [grep {($_->description || '') !~ /$bam_file/} @{$tr->get_all_Attributes}];
    $tr->edits_enabled($edits_enabled_bak);

    return;
  }

  return $tr;
}



## TREE METHODS
###############


=head2 transcript_tree

  Example    : $tree = $as->transcript_tree();
  Description: Gets transcript tree for this annotation source.
  Returntype : Bio::EnsEMBL::VEP::TranscriptTree
  Exceptions : none
  Caller     : get_nearest()
  Status     : Stable

=cut

sub transcript_tree {
  my $self = shift;

  unless(exists($self->{transcript_tree})) {
    $self->{transcript_tree} = Bio::EnsEMBL::VEP::TranscriptTree->new({
      config => $self->config,
      annotation_source => $self
    });
  }

  return $self->{transcript_tree};
}


=head2 populate_tree

  Arg 1      : Bio::EnsEMBL::VEP::TranscriptTree $tree
  Example    : $as->populate_tree($tree);
  Description: Populates transcript tree with data from this annotation source.
  Returntype : none
  Exceptions : none
  Caller     : TranscriptTree::new()
  Status     : Stable

=cut

sub populate_tree {
  my ($self, $tree) = @_;

  # insert into tree from file
  open TR, $self->tree_file or throw("ERROR: Could not read from tree file: $!");
  $self->_tree_insert_file_line($tree, $_) while (<TR>);
  close TR;
}


=head2 get_nearest

  Arg 1      : Bio::EnsEMBL::Variation::BaseVariationFeature $vf
  Arg 2      : string $type (transcript, gene or symbol)
  Example    : $nearest_gene = $as->get_nearest($vf, 'symbol');
  Description: Gets ID (one of transcript, gene, symbol) of the nearest
               transcript to the given variant.
  Returntype : string
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub get_nearest {
  my ($self, $vf, $type) = @_;

  throw("ERROR: No type supplied for --nearest\n") unless $type;

  my @return;
  my %seen;
  foreach my $nearest(@{$self->transcript_tree->nearest($vf->{chr}, $vf->{start}, $vf->{end})}) {
    throw("ERROR: Invalid type \"$type\" for --nearest\n") unless exists($nearest->{$type});
    next unless my $value = $nearest->{$type};
    push @return, $value unless $seen{$value};
    $seen{$value} = 1;
  }

  return \@return;
}


=head2 _tree_coords_filename

  Example    : $file = $as->_tree_coords_filename();
  Description: Get filename to read/write transcript coordinates used to
               populate transcript tree.
  Returntype : string
  Exceptions : none
  Caller     : populate_tree()
  Status     : Stable

=cut

sub _tree_coords_filename {
  return $_[0]->dir.'/transcript_gene_tss.txt';
}


=head2 _tree_file_data

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Example    : $data = $as->_tree_file_data($tr);
  Description: Gets arrayref of data for writing to tree file from a transcript.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : tree_file()
  Status     : Stable

=cut

sub _tree_file_data {
  my ($self, $tr) = @_;
  my $tss = $tr->seq_region_strand == 1 ? $tr->seq_region_start : $tr->seq_region_end;
  return [$tss, $tr->stable_id, $tr->{_gene_stable_id}, ($tr->{_gene_symbol} || $tr->{_gene_stable_id})];
}


=head2 _tree_insert_file_line

  Arg 1      : Bio::EnsEMBL::VEP::TranscriptTree
  Example    : $as->_tree_insert_file_line($tree, $line);
  Description: Insert a line of data read from tree file to the transcript tree.
  Returntype : none
  Exceptions : none
  Caller     : populate_tree()
  Status     : Stable

=cut

sub _tree_insert_file_line {
  my ($self, $tree, $line) = @_;
  chomp($line);
  my ($c, $tss, $tr_stable_id, $gene_stable_id, $gene_symbol) = split("\t", $line);
  $tree->insert(
    $c, $tss, $tss,
    {
      s          => $tss,
      e          => $tss,
      transcript => $tr_stable_id,
      gene       => $gene_stable_id,
      symbol     => $gene_symbol
    }
  );
}

1;
