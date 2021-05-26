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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF - parent class for GFF and GTF annotation sources

=head1 SYNOPSIS

Should not be invoked directly

=head1 DESCRIPTION

Base class for GFF and GTF custom file annotation sources.

=head1 METHODS

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
use Bio::EnsEMBL::Attribute;

use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
$Data::Dumper::Indent = 1;

use base qw(
  Bio::EnsEMBL::VEP::AnnotationType::Transcript
  Bio::EnsEMBL::VEP::AnnotationSource::File
);


=head2 new

  Arg 1      : hashref $args
               {
                 config            => Bio::EnsEMBL::VEP::Config $config,
                 file              => string $filename,
                 type              => 'exact',
                 short_name        => (optional) string $short_name,
                 transcript_filter => (optional) string $filter,
                 bam               => (optional) string $bam_file
               }
  Example    : Not invoked directly
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::File::BaseGXF
  Exceptions : throws if:
                - no database or FASTA available to read sequence
                - type is not "exact"
               warns if --nearest specified, not currently available
  Caller     : Bio::EnsEMBL::VEP::AnnotationSource::File->new()
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->Bio::EnsEMBL::VEP::AnnotationSource::File::new(@_);

  # requires sequence
  throw("ERROR: GXF annotation requires either database access (--database or --cache) or a FASTA file (--fasta)")
    unless $self->param('fasta') or ($self->param('cache') && !$self->param('offline')) or $self->param('database');

  # no point doing exact
  throw("ERROR: GXF annotation sources cannot be used as \"exact\" custom annotation type\n") if $self->type eq 'exact';

  $self->add_shortcuts(['use_transcript_ref']);

  if($self->param('nearest')) {
    $self->warning_msg("WARNING: --nearest is not currently compatible with GXF annotation sources");
  }

  $self->{cache_region_size} = 1e6;

  return $self;
}


=head2 get_features_by_regions_uncached

  Arg 1      : arrayref $regions
  Example    : $transcripts = $as->get_features_by_regions_uncached($regions)
  Description: Gets all transcripts overlapping the given set of regions. See
               Bio::EnsEMBL::VEP::AnnotationSource::get_all_regions_by_InputBuffer()
               for information about regions. Transcripts are not yet EnsEMBL
               objects; they are converted if lazy_load_transcript() is called.
               This means only transcripts that actually overlap input variants
               undergo the (time-costly) object conversion process.
  Returntype : arrayref of transcript record hashes
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

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
      ($s * $cache_region_size) + 1,
      ($s + 1) * $cache_region_size
    );

    $cache->{$c}->{$s} = $features;

    push @return, @$features;
  }

  return \@return;
}


=head2 _get_transcripts_by_coords

  Arg 1      : string $chr
  Arg 2      : int $start
  Arg 3      : int $end
  Example    : $transcripts = $as->_get_transcripts_by_coords($chr, $start, $end)
  Description: Gets all transcripts overlapping the given chromosomal region.
  Returntype : arrayref of of transcript record hashes
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub _get_transcripts_by_coords {
  my $self = shift;
  return $self->_create_transcripts($self->_get_records_by_coords(@_));
}


=head2 _get_records_by_coords

  Arg 1      : string $chr
  Arg 2      : int $start
  Arg 3      : int $end
  Arg 4      : (optional) int $no_rescan
  Example    : $records = $as->_get_records_by_coords($chr, $start, $end)
  Description: Gets all record hashrefs from file overlapping given coordinates.
               Runs iteratively if a record overhangs either end to pick up
               sub-records (e.g. exons). Set $no_rescan to prevent iterative
               running:
                - $no_rescan = 1 : don't expand search to lower coord values
                - $no_rescan = 2 : don't expand search to higher coord values
  Returntype : arrayref of of transcript record hashes
  Exceptions : none
  Caller     : _get_transcripts_by_coords()
  Status     : Stable

=cut

sub _get_records_by_coords {
  my ($self, $c, $s, $e, $no_rescan) = @_;

  my $parser = $self->parser();
  if ($parser->seek($c, $s - 1, $e + 1)) {
    $parser->next();
  } else {
    return;
  }

  my $include = $self->include_feature_types;

  my ($min, $max) = (1e10, 0);
  my @records;
  
  while($parser->{record} && $parser->get_start <= $e) {

    if($include->{$parser->get_type}) {
      my ($r_start, $r_end) = ($parser->get_start, $parser->get_end);
      
      # if this is a rescan, only change min/max in the appropriate "direction"
      if($no_rescan) {
        $min = $r_start if $r_start < $min && $no_rescan == 1;
        $max = $r_end if $r_end > $max && $no_rescan == 2;
      }
      else {
        $min = $r_start if $r_start < $min;
        $max = $r_end if $r_end > $max;
      }

      push @records, $self->_record_to_hash();
    }

    $parser->next();
  }

  # we have to be a bit clever, because sub-features e.g. exons may not come back
  # as they dont overlap our input coords, even if the parent feature does!
  # so rerun this method but set no_rescan to an int value representing the "direction"
  # so we don't keep re-iterating up/down/up/down
  if($min < $s) {
    unshift @records, @{$self->_get_records_by_coords($c, $min, $s - 1, 1)};
  }
  if($max > $e) {
    push @records, @{$self->_get_records_by_coords($c, $e + 1, $max, 2)};
  }

  return \@records;
}


=head2 _record_to_hash

  Example    : $record_hash = $as->_record_to_hash()
  Description: Gets a hashref representing the current parser record; contains
               all pertinent record information and an md5 key to ease unique
               sorting.
  Returntype : hashref
  Exceptions : none
  Caller     : _get_records_by_coords()
  Status     : Stable

=cut

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


=head2 _create_transcripts

  Arg 1      : arrayref $record_hashes
  Example    : $transcripts = $as->_create_transcripts($record_hashes)
  Description: Takes a list of record hashes as read from the source file
               and converts them to transcript objects. This is done by
               creating a "tree" representing the parent/child structure
               of the records in the GFF/GTF
               (typically gene->transcript->exon), then compiling these
               into EnsEMBL objects (Gene, Transcript, Exon).
  Returntype : arrayref of transcript record hashes
  Exceptions : none
  Caller     : _get_records_by_coords()
  Status     : Stable

=cut

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
      $tr_record->{_gene_record} = $gene;
      push @transcripts, $tr_record;
    }

    # avoid a circular ref
    delete($gene->{_children}) if $gene;
  }

  return \@transcripts;
}


=head2 _get_parent_child_structure

  Arg 1      : arrayref $record_hashes
  Example    : $tree = $as->_get_parent_child_structure($record_hashes)
  Description: Takes a list of record hashes as read from the source file
               and converts them to a "tree" representing the
               parent/child structure of the records in the GFF/GTF
               (typically gene->transcript->exon)
  Returntype : hashref
  Exceptions : throws if record found that references non-existent parent
  Caller     : _create_transcripts()
  Status     : Stable

=cut

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
  my %parent_children;
  foreach my $record(@$records) {
    my $record_id = $self->_record_get_id($record);

    foreach my $parent_id (@{$self->_record_get_parent_ids($record)}) {
      next unless ($parent_id);
      if(my $top_level = $top_level{$parent_id}) {
        push @{$top_level->{_children}}, $record;
        # record id can be shared across multiple CDS lines, so key on position as well as id
        $parent_children{$parent_id}{$record_id}{$record->{start}} = 1;
      }
      else {
        push @orphans, $record;
      }
    }

    $top_level{$record_id} = $record;
  }

  # now deal with orphans
  my %parents_not_found = ();

  foreach my $record(@orphans) {
    my $record_id = $self->_record_get_id($record);

    foreach my $parent_id (@{$self->_record_get_parent_ids($record)}) {
      # Check if the relation parent/record has already been found: use case where a record is related to more than one parent record
      # e.g. an exon linked to several transcripts
      next if ($parent_children{$parent_id}{$record_id}{$record->{start}});

      if(my $top_level = $top_level{$parent_id}) {
        push @{$top_level->{_children}}, $record;
        $parent_children{$parent_id}{$record_id}{$record->{start}} = 1;
      }
      else {
        $parents_not_found{$parent_id} = 1;
      }
    }
  }

  $self->warning_msg(
    "WARNING: Parent entries with the following IDs were not found or skipped due to invalid types: ".
    join(", ", keys %parents_not_found)
  ) if scalar keys %parents_not_found;

  # now prune the structure so we're left with only true top level records
  foreach my $entry_id (keys(%top_level)) {
    my $parent_ids = $self->_record_get_parent_ids($top_level{$entry_id});
    if (scalar(@$parent_ids) && $parent_ids->[0] ne '') {
      delete $top_level{$entry_id};
    }
  }
  return \%top_level;
}


=head2 _create_transcript

  Arg 1      : hashref $transcript_record_hash
  Arg 2      : hashref $gene_record_hash
  Example    : $transcript = $as->_create_transcript($tr_hash, $gene_hash)
  Description: Takes a structured transcript record hash and creates a
               Bio::EnsEMBL::Transcript object from it.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throws if record has child records not of type "exon", "cds" or "stop_codon"
               warns if:
                - unable to determine biotype
                - unable to add Exon object (typically if coordinates overlap)
  Caller     : lazy_load_transcript()
  Status     : Stable

=cut

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
    $self->warning_msg("WARNING: Unable to determine biotype of $id");
    return;  
  }

  my ($tr_start, $tr_end, $tr_strand) = (
    $tr_record->{start},
    $tr_record->{end},
    $tr_record->{strand}
  );

  my $tr = Bio::EnsEMBL::Transcript->new_fast({
    stable_id => $id,
    biotype   => $biotype,
    slice     => $slice,
    strand    => $tr_strand,
    version   => 1,
    dbID      => $self->{_tr_dbID}++,
    start     => $tr_start,
    end       => $tr_end,
  });

  # making a call to seq here means the sequence is cached in memory
  # this means we don't need to repeatedly fetch for each exon
  # the caching is done elsewhere so we don't need to worry any further here
  $self->_quick_sub_Slice($slice, $tr_start, $tr_end, $tr_strand)->seq;

  $self->_add_identifiers($tr, $tr_record, $gene_record);

  # separate exons and cds entries
  my (%cdss, @exons);

  foreach my $child(@{$tr_record->{_children}}) {
    my $type = lc($child->{type});

    if($type eq 'exon') {
      push @exons, $child;
    }
    elsif($type eq 'cds' || $type eq 'stop_codon') {
      $cdss{$child->{start}} = $child;
    }
    else {
      throw("ERROR: Transcript has unexpected type of child record: ".Dumper($child)."\n");
    }
  }

  # check for exons for protein_coding biotype
  if ($biotype eq 'protein_coding' && scalar @exons == 0){
    $self->warning_msg("WARNING: No exons found for protein_coding transcript $id");
    return;
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
    my ($s, $e, $str) = ($exon_record->{start}, $exon_record->{end}, $exon_record->{strand});

    my $phase = -1;
    my @cds_records;

    if(my @cds_starts = grep {overlap($s, $e, $_, $_)} keys %cdss) {
      if(@cds_starts > 1) {
        if($tr_record->{strand} > 0) {
          @cds_starts = sort {$a <=> $b} @cds_starts;
        }
        else {
          @cds_starts = sort {$b <=> $a} @cds_starts;
        }
      }

      @cds_records = map {delete $cdss{$_}} @cds_starts;
      push @ordered_cdss, @cds_records;
      $phase = $self->_convert_phase($cds_records[0]->{phase});
    }

    my $exon = Bio::EnsEMBL::Exon->new(
      -START  => $s,
      -END    => $e,
      -STRAND => $str,
      -SLICE  => $slice,
      -PHASE  => $phase,
    );

    $exon->{_seq_cache} = $self->_quick_sub_Slice($slice, $s, $e, $str)->seq;

    # log a pointer to the exon on the cds record
    $_->{_exon} = $exon for @cds_records;

    # add it to the transcript
    # sometimes this can fail if the coordinates overlap
    eval {$tr->add_Exon($exon, $exon_record->{attributes}->{rank});};
    if($@) {
      $self->warning_msg("WARNING: Failed to add exon to transcript ".$tr->stable_id."\n$@");
      return;
    }
  }

  $self->_add_translation($tr, \@ordered_cdss) if @ordered_cdss;

  return $tr;
}


=head2 _quick_sub_Slice

  Arg 1      : Bio::EnsEMBL::Slice
  Arg 2      : int $start
  Arg 3      : int $end
  Arg 4      : int $strand
  Example    : $sub = $as->_quick_sub_Slice($slice, $start, $end, $strand)
  Description: Faster way of creating sub slices that $slice->sub_Slice()
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : _create_transcript()
  Status     : Stable

=cut

sub _quick_sub_Slice {
  my ($self, $slice, $start, $end, $strand) = @_;

  my %new = %{$slice};
  $new{start}  = $start;
  $new{end}    = $end;
  $new{strand} = $strand;

  bless(\%new, ref($slice));

  return \%new;
}


=head2 _add_identifiers

  Arg 1      : Bio::EnsEMBL::Transcript
  Arg 2      : hashref $transcript_record_hash
  Arg 3      : hashref $gene_record_hash
  Example    : $as->_add_identifiers($tr, $tr_record, $gene_record)
  Description: Add any identifiers we can find in the transcript and/or gene
               record to the transcript object. Currently this is only the
               gene stable_id and gene symbol.
  Returntype : none
  Exceptions : none
  Caller     : _create_transcript()
  Status     : Stable

=cut

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

  # add TSL
  if(my $tsl = $tr_record->{attributes}->{transcript_support_level} || $tr_record->{attributes}->{tsl}) {
    push @{$tr->{attributes}}, Bio::EnsEMBL::Attribute->new_fast({
      code => 'TSL',
      value => 'tsl'.$tsl,
    });
  }
}


=head2 _add_translation

  Arg 1      : Bio::EnsEMBL::Transcript
  Arg 2      : arrayref $cds_record_hashes
  Example    : $translation = $as->_add_translation($tr, $tr_record, $gene_record)
  Description: Create the translation object for this transcript.
  Returntype : Bio::EnsEMBL::Translation
  Exceptions : none
  Caller     : _create_transcript()
  Status     : Stable

=cut

sub _add_translation {
  my ($self, $tr, $ordered_cdss) = @_;

  # get ID if available
  my $protein_id = $ordered_cdss->[0]->{attributes}->{protein_id} || undef;

  # copy protein ID to transcript object for outputfactory
  $tr->{_protein} = $protein_id if $protein_id;

  # create translation object
  my $translation = $tr->{translation} ||= Bio::EnsEMBL::Translation->new(
    -TRANSCRIPT => $tr,
    -VERSION    => 1,
    -STABLE_ID  => $protein_id,
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

  my $chr = $tr->{slice}->seq_region_name;

  # translate
  # we have to delete slice otherwise the API tries to look up codon tables etc
  # from a non-existent adaptor
  my $slice = delete($tr->{slice});
  $tr->{_variation_effect_feature_cache}->{peptide} = $translation->seq;
  $tr->{_variation_effect_feature_cache}->{codon_table} = 1;
  if($chr eq 'MT') {
    $tr->{_variation_effect_feature_cache}->{codon_table} = 2;
  }

  $tr->{slice} = $slice;

  return $translation;
}


=head2 _record_is_gene

  Arg 1      : hashref $record_hash
  Example    : $is_gene = $as->_record_is_gene($record_hash)
  Description: Establish if a given record is a gene-level record.
  Returntype : bool
  Exceptions : none
  Caller     : _create_transcript()
  Status     : Stable

=cut

sub _record_is_gene {
  my ($self, $record) = @_;

  return 1 if $record->{type} =~ /gene/;
  return 1 if $record->{type} eq 'processed_transcript';
  return 1 if $record->{type} eq 'RNA';

  return 0;
}


=head2 _convert_phase

  Arg 1      : int $phase
  Example    : $ensembl_phase = $as->_convert_phase($gff_phase)
  Description: Converts GFF exon phase to Ensembl phase.
               0 => 0, 1 => 2, 2 => 1
  Returntype : bool
  Exceptions : none
  Caller     : _create_transcript()
  Status     : Stable

=cut

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


=head2 lazy_load_transcript

  Arg 1      : hashref $transcript_record_hash
  Example    : $tr = $as->lazy_load_transcript($tr_record_hash)
  Description: Converts given structured transcript record hash to Ensembl
               object. This will only be called when the object is about to
               be used, rather than when it is read from disk, saving time.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub lazy_load_transcript {
  my $self = shift;
  my $feature = shift;

  $feature->{object} = $self->_create_transcript($feature, $feature->{_gene_record}) if not exists($feature->{object});

  return $feature->{object};
}


=head2 merge_features

  Arg 1      : arrayref $transcript_record_hashes
  Example    : $unique_trs = $as->merge_features($tr_record_hashes)
  Description: Return a unique list of transcript record hashes. Unique
               sorting done on md5 key added when record was created.
  Returntype : arrayref
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub merge_features {
  my ($self, $features) = @_;
  my (@return, %seen);

  foreach my $f(@$features) {
    unless($seen{$f->{md5}}) {
      push @return, $f;
      $seen{$f->{md5}} = 1;  
    }
  }

  return \@return;
}

1;
