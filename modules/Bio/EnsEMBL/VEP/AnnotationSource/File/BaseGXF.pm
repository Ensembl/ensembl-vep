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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF - parent class for GFF and GTF annotation sources

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
$Data::Dumper::Indent = 1;

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::BaseTranscript
  Bio::EnsEMBL::VEP::AnnotationSource::File
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # requires sequence
  throw("ERROR: GXF annotation requires either database access (--database or --cache) or a FASTA file (--fasta)")
    unless $self->param('fasta') or $self->param('cache') or $self->param('database');

  $self->{cache_region_size} = 1e6;

  return $self;
}

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;
  my $cache_region_size = $self->{cache_region_size};

  my @return;

  foreach my $region(@{$regions}) {
    my ($c, $s) = @$region;

    my $features = $self->_get_transcripts_by_coords(
      $c,
      $s * $cache_region_size,
      (($s + 1) * $cache_region_size) - 1
    );

    $cache->{$c}->{$s} = $features;

    push @return, @$features;
  }

  return \@return;
}

sub _get_transcripts_by_coords {
  my $self = shift;
  return $self->_create_transcripts($self->_get_records_by_coords(@_));
}

sub _get_records_by_coords {
  my ($self, $c, $s, $e, $no_rescan) = @_;

  my $parser = $self->parser();
  $parser->seek($c, $s - 1, $e + 1);
  $parser->next();

  my $include = $self->include_feature_types;

  my ($min, $max) = (1e10, 0);
  my @records;
  
  while($parser->{record} && $parser->get_start <= $e) {

    if($include->{$parser->get_type}) {
      my ($r_start, $r_end) = ($parser->get_start, $parser->get_end);
      
      # if this is a rescan, don't change min/max so we don't get repeated iterations
      unless($no_rescan) {
        $min = $r_start if $r_start < $min;
        $max = $r_end if $r_end > $max;
      }

      push @records, $self->_record_to_hash();
    }

    $parser->next();
  }

  # we have to be a bit clever, because sub-features e.g. exons may not come back
  # as they dont overlap our input coords, even if the parent feature does!
  # so rerun this method but set no_rescan so we don't keep re-iterating
  if($min < $s) {
    unshift @records, @{$self->_get_records_by_coords($c, $min, $s - 1, 1)};
  }
  if($max > $e) {
    push @records, @{$self->_get_records_by_coords($c, $e + 1, $max, 1)};
  }

  return \@records;
}

sub _record_to_hash {
  my $self = shift;
  my $parser = $self->parser;

  return {
    md5        => md5_hex($parser->{current_block}),
    chr        => $parser->get_seqname,
    start      => $parser->get_start,
    end        => $parser->get_end,
    strand     => $parser->get_strand,
    phase      => $parser->get_phase,
    source     => $parser->get_source,
    type       => $parser->get_type,
    attributes => $parser->get_attributes,
  };
}

sub _create_transcripts {
  my $self = shift;
  my $records = shift;

  # we want a structure relating sub-features to their parents
  my $top_level = $self->_get_parent_child_structure($records);

  my @transcripts;

  foreach my $record(values %$top_level) {

    my $gene;
    my @gene_transcripts;

    if($self->_record_is_gene($record)) {
      $gene = $record;

      # some RefSeq gene records have no children...
      @gene_transcripts = @{$record->{_children} || []};
    }
    else {
      @gene_transcripts = ($record);
    }

    foreach my $tr_record(@gene_transcripts) {
      my $tr = $self->_create_transcript($tr_record, $gene);
      push @transcripts, $tr if $tr;
    }
  }

  return \@transcripts;
}

sub _get_parent_child_structure {
  my $self = shift;
  my $records = shift;

  # unique sort using md5
  my @new = ();
  my %seen = ();

  for(@$records) {
    next if $seen{$_->{md5}};
    push @new, $_;
    $seen{$_->{md5}} = 1;
  }

  $records = \@new;

  my %top_level;

  # sub-features frequently won't appear in order after their parents
  # so we need to track orphans and add them later
  my @orphans;

  foreach my $record(@$records) {
    if(my $parent_id = $self->_record_get_parent_id($record)) {
      if(my $top_level = $top_level{$parent_id}) {
        push @{$top_level->{_children}}, $record;
      }
      else {
        push @orphans, $record;
      }
    }

    $top_level{$self->_record_get_id($record)} = $record;
  }

  # now deal with orphans
  foreach my $record(@orphans) {
    my $parent_id = $self->_record_get_parent_id($record);

    if(my $top_level = $top_level{$parent_id}) {
      push @{$top_level->{_children}}, $record;
    }
    else {
      throw("ERROR: Parent record for the following not found:\n".Dumper($record)."\n");
    }
  }

  # now prune the structure so we're left with only true top level records
  delete $top_level{$_} for grep {$self->_record_get_parent_id($top_level{$_})} keys %top_level;

  return \%top_level;
}

sub _create_transcript {
  my $self = shift;
  my $tr_record = shift;
  my $gene_record = shift;

  return unless $tr_record->{_children};

  my $id = $tr_record->{attributes}->{transcript_id} || $self->_record_get_id($tr_record);
  $id =~ s/^(gene|transcript)://i;

  my $slice = $self->get_slice($tr_record->{chr});

  my $biotype = $self->_record_get_biotype($tr_record, $gene_record);
  unless($biotype) {
    warn("WARNING: Unable to determine biotype of $id");
    return;  
  }  

  my $tr = Bio::EnsEMBL::Transcript->new(
    -STABLE_ID => $id,
    -BIOTYPE   => $biotype,
    -SLICE     => $slice,
    -STRAND    => $tr_record->{strand},
    -VERSION   => 1,
    -dbID      => $self->{_tr_dbID}++,
  );

  $self->_add_identifiers($tr, $tr_record, $gene_record);

  # separate exons and cds entries
  my (@cdss, @exons);
  foreach my $child(@{$tr_record->{_children}}) {
    my $type = lc($child->{type});

    if($type eq 'exon') {
      push @exons, $child;
    }
    elsif($type eq 'cds') {
      push @cdss, $child;
    }
    else {
      throw("ERROR: Transcript has unexpected type of child record: ".Dumper($child)."\n");
    }
  }

  # sort exons
  if($tr_record->{strand} > 0) {
    @exons = sort {$a->{start} <=> $b->{start}} @exons;
  }
  else {
    @exons = sort {$b->{start} <=> $a->{start}} @exons;
  }

  # now create exon objects and add them to the transcript
  my @ordered_cdss;

  foreach my $exon_record(@exons) {
    my ($s, $e) = ($exon_record->{start}, $exon_record->{end});

    my $phase = -1;
    my $cds_record;
    if(($cds_record) = grep {overlap($s, $e, $_->{start}, $_->{end})} @cdss) {
      push @ordered_cdss, $cds_record;
      $phase = $self->_convert_phase($cds_record->{phase});
    }

    my $exon = Bio::EnsEMBL::Exon->new(
      -START  => $s,
      -END    => $e,
      -STRAND => $exon_record->{strand},
      -SLICE  => $slice,
      -PHASE  => $phase,
    );

    $exon->{_seq_cache} = $exon->feature_Slice->seq;

    # log a pointer to the exon on the cds record
    $cds_record->{_exon} = $exon if $cds_record;

    # add it to the transcript
    # sometimes this can fail if the coordinates overlap
    eval {$tr->add_Exon($exon);};
    if($@) {
      warn("WARNING: Failed to add exon to transcript ".$tr->stable_id."\n$@");
      return;
    }
  }

  $self->_add_translation($tr, \@ordered_cdss) if @ordered_cdss;

  return $tr;
}

sub _add_identifiers {
  my ($self, $tr, $tr_record, $gene_record) = @_;

  $tr->{_source_cache} = $self->short_name;

  # get gene ID
  $tr->{_gene_stable_id} = $tr_record->{attributes}->{gene_id};

  if(!$tr->{_gene_stable_id}) {
    foreach my $pair(split(',', $tr_record->{attributes}->{dbxref} || $tr_record->{attributes}->{Dbxref} || '')) {
      my ($k, $v) = split(':', $pair);
      if($k eq 'GeneID') {
        $tr->{_gene_stable_id} = $v;
        last;
      }
    }
  }

  if(!$tr->{_gene_stable_id} && $gene_record) {
    if($gene_record->{attributes}->{gene_id}) {
      $tr->{_gene_stable_id} = $gene_record->{attributes}->{gene_id};
    }
    else {
      $tr->{_gene_stable_id} = $self->_record_get_id($gene_record);
      $tr->{_gene_stable_id} =~ s/^(gene|transcript)\://i;
    }
  }

  # try and get the gene symbol
  for my $key(qw(gene_name gene)) {
    if($tr_record->{attributes}->{$key}) {
      $tr->{_gene_symbol} = $tr_record->{attributes}->{$key};
      last;
    }
  }

  if(!$tr->{_gene_symbol} && $gene_record) {
    $tr->{_gene_symbol} = $gene_record->{attributes}->{name} || $gene_record->{attributes}->{Name};
  }
}

sub _add_translation {
  my ($self, $tr, $ordered_cdss) = @_;

  # create translation object
  my $translation = $tr->{translation} ||= Bio::EnsEMBL::Translation->new(
    -TRANSCRIPT => $tr,
    -VERSION    => 1,
  );
  $translation->{transcript} = $tr;
  weaken($translation->{transcript});

  # we need to use the first and last cds record to set the translation start and ends
  my $offset;

  if($ordered_cdss->[0]->{strand} > 0) {
    $offset = ($ordered_cdss->[0]->{start} - $ordered_cdss->[0]->{_exon}->start) + 1;
  }
  else {
    $offset = ($ordered_cdss->[0]->{_exon}->end - $ordered_cdss->[0]->{end}) + 1;
  }
  $translation->start($offset);
  $translation->start_Exon($ordered_cdss->[0]->{_exon});

  if($ordered_cdss->[-1]->{strand} > 0) {
    $offset = ($ordered_cdss->[-1]->{end} - $ordered_cdss->[-1]->{_exon}->start) + 1;
  }
  else {
    $offset = ($ordered_cdss->[-1]->{_exon}->end - $ordered_cdss->[-1]->{start}) + 1;
  }
  $translation->end($offset);
  $translation->end_Exon($ordered_cdss->[-1]->{_exon});

  # translate
  # we have to delete slice otherwise the API tries to look up codon tables etc
  # from a non-existent adaptor
  my $slice = delete($tr->{slice});
  $tr->{_variation_effect_feature_cache}->{peptide} = $translation->seq;
  $tr->{_variation_effect_feature_cache}->{codon_table} = 1;
  $tr->{slice} = $slice;

  return $translation;
}

sub _record_get_parent_id {
  my ($self, $record) = @_;

  if(!exists($record->{_parent_id})) {
    my $attributes = $record->{attributes};
    $record->{_parent_id} = $attributes->{Parent} || $attributes->{parent};
  }

  return $record->{_parent_id};
}

sub _record_get_id {
  my ($self, $record) = @_;

  if(!exists($record->{_id})) {
    my $attributes = $record->{attributes};
    $record->{_id} = $attributes->{ID} || $attributes->{Name} || $attributes->{id} || $attributes->{name};
  }

  return $record->{_id};
}

sub _record_get_biotype {
  my ($self, $record, $gene_record) = @_;

  if(!exists($record->{_biotype})) {

    # Ensembl-y GFFs have biotype as an attribute
    my $biotype = $record->{attributes}->{biotype};

    # others we need to (guess) work it out
    if(!$biotype) {
      my $type = lc($record->{type});

      if($type eq 'mrna') {
        $biotype = 'protein_coding';
      }
      elsif($type eq 'ncrna') {
        $biotype = $record->{attributes}->{ncrna_class};
      }
      elsif($type =~ /^([a-z]+)_gene_segment$/) {
        $biotype = 'IG_'.uc($1).'_gene';
      }
      elsif($gene_record && ($gene_record->{attributes}->{description} || '') =~ /^microRNA/) {
        $biotype = 'miRNA';
      }
      elsif($record->{attributes}->{gbkey}) {
        $biotype = $record->{attributes}->{gbkey};
      }
    }

    $record->{_biotype} = $biotype;
  }

  return $record->{_biotype};
}

sub _record_is_gene {
  my ($self, $record) = @_;

  return 1 if $record->{type} =~ /gene/;
  return 1 if $record->{type} eq 'processed_transcript';
  return 1 if $record->{type} eq 'RNA';

  return 0;
}

# converts GTF/GFF phase to Ensembl phase
# basically 0 => 0, 1 => 2, 2 => 1; why? Because.
sub _convert_phase {
  my ($self, $phase) = @_;

  if($phase == 1) {
    $phase = 2;
  }
  elsif($phase == 2) {
    $phase = 1;
  }

  return $phase;
}

1;