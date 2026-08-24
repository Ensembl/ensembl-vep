=head1 LICENSE

Copyright [2016-2026] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat - GFF-based regulatory feature annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
  config => $config,
  file   => "regulatory_features.gff3.gz"
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Reads regulatory features from a tabix-indexed GFF3 file and produces genuine
regulatory consequences, as an alternative to the funcgen database
(Bio::EnsEMBL::VEP::AnnotationSource::Database::RegFeat) or the VEP cache
(Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat).

We must create RegulatoryFeature objects because regulatory consequence
assignment only works for that class. It is gated on the object's Perl class:
_get_oc_list blesses a dummy into the OverlapConsequence's feature_class and
calls isa() on it, so a generic Bio::EnsEMBL::Feature yields no regulatory
consequences and no error (see
Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele::_get_oc_list).

Only eight fields survive into the VEP cache (see
Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Regulation::clean_regfeat), so the
objects built here carry exactly those and nothing more.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat;

use Scalar::Util qw(weaken);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::MotifFeature;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;
use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::IO::Parser::GFF3Tabix;

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::File
  Bio::EnsEMBL::VEP::AnnotationType::RegFeat
);

# GFF column 3 values accepted as regulatory features. These are all valid SO
# terms and are used verbatim as the feature_type, which becomes BIOTYPE in
# output - the cache stores a bare SO string here too, so this is exact parity.
our %INCLUDE_FEATURE_TYPES = map {$_ => 1} qw(
  promoter
  enhancer
  CTCF_binding_site
  open_chromatin_region
  promoter_flanking_region
  TF_binding_site
  regulatory_region
  EMAR
);

# only promoters carry extended_start/extended_end
our %HAS_EXTENDED_BOUNDS = map {$_ => 1} qw(
  promoter
);

# margin to pad the query when extended_promoters is active, so a promoter whose
# core lies in a neighbouring region but whose extended bounds reach into this
# one is still found.
our $EXTENDED_PROMOTER_MARGIN = 10_000;

# GFF column 3 values accepted when this source is reading a motif GFF
# (--regulatory_gff motifs=). Kept separate from the regulatory set so a motif
# file and a regulatory file each accept only their own records.
our %MOTIF_INCLUDE_FEATURE_TYPES = map {$_ => 1} qw(
  TF_binding_site
);

# Length given to a fabricated slice when no real one is available. Offline
# there is no way to know the true sequence length, so this only has to be large
# enough never to truncate a feature; it exceeds the longest known chromosome.
our $FAKE_SLICE_LENGTH = 2_000_000_000;


=head2 new

  Arg 1      : hashref $args
               {
                 config             => Bio::EnsEMBL::VEP::Config $config,
                 file               => string $filename,
                 extended_promoters => (optional) bool,
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat
  Exceptions : throws if no file is given
  Caller     : AnnotationSourceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  # bypass AnnotationSource::File::new, which dispatches on {format} to a
  # format subclass; this class is selected directly by --regulatory_gff
  my $self = $class->Bio::EnsEMBL::VEP::AnnotationSource::new(@_);

  my $hashref = $_[0];

  throw("ERROR: No file given\n") unless $hashref->{file};
  $self->file($hashref->{file});

  # A motif source (--regulatory_gff motifs=) builds MotifFeatures; otherwise the
  # source builds RegulatoryFeatures. Everything else - parsing, dedup, slices,
  # identifier handling - is shared, mirroring how Cache::RegFeat and
  # Database::RegFeat each emit both feature classes from one class.
  $self->{is_motif} = $hashref->{motif} ? 1 : 0;

  $self->short_name($hashref->{short_name} ||
    ($self->{is_motif} ? 'MotifFeatures' : 'RegulatoryFeatures'));
  $self->type('overlap');

  # regulatory-only options; a motif GFF carries no extended bounds
  unless($self->{is_motif}) {
    $self->add_shortcuts([qw(extended_promoters)]);

    # A GFF carries no epigenome activity, so --cell_type cannot be honoured.
    # Warn and ignore rather than fail: CELL_TYPE is simply empty for
    # GFF-derived features, and the warning makes that explicit.
    #
    # warning_msg deduplicates per object, not per run, so a second regulatory
    # instance (EMARs) would repeat it. The adaptor sets quiet_cell_type on
    # secondary instances so the warning is emitted once.
    if(my $ct = $self->param('cell_type')) {
      $self->warning_msg(
        "WARNING: --cell_type is ignored for regulatory features from ".
        "--regulatory_gff; a GFF carries no epigenome activity data"
      ) if (ref($ct) eq 'ARRAY' ? scalar @$ct : $ct) && !$hashref->{quiet_cell_type};
    }
  }

  $self->{cache_region_size} = 1e6;

  return $self;
}


=head2 parser

  Example    : $parser = $as->parser();
  Description: Get ensembl-io parser to read from file.
  Returntype : Bio::EnsEMBL::IO::Parser::GFF3Tabix
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::GFF3Tabix->open(
    $self->file, must_parse_metadata => 0
  );
}


=head2 include_feature_types

  Example    : $types = $as->include_feature_types();
  Description: Get hashref of GFF column 3 values this class accepts.
  Returntype : hashref
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub include_feature_types {
  my $self = shift;
  return $self->{is_motif}
    ? \%MOTIF_INCLUDE_FEATURE_TYPES
    : \%INCLUDE_FEATURE_TYPES;
}


=head2 get_features_by_regions_uncached

  Arg 1      : arrayref $regions
  Example    : $features = $as->get_features_by_regions_uncached($regions)
  Description: Gets all regulatory features overlapping the given set of regions.
               See Bio::EnsEMBL::VEP::AnnotationSource::get_all_regions_by_InputBuffer()
               for information about regions.
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::RegulatoryFeature
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

    my $features = $self->_get_regfeats_by_coords(
      $c,
      ($s * $cache_region_size) + 1,
      ($s + 1) * $cache_region_size
    );

    $cache->{$c}->{$s} = $features;

    push @return, @$features;
  }

  return \@return;
}


=head2 _get_regfeats_by_coords

  Arg 1      : string $chr
  Arg 2      : int $start
  Arg 3      : int $end
  Example    : $features = $as->_get_regfeats_by_coords($chr, $start, $end)
  Description: Gets all regulatory features overlapping the given region.
               Unlike BaseGXF this does no parent/child rescanning: regulatory
               features are flat intervals with no sub-features. The query is
               padded by $EXTENDED_PROMOTER_MARGIN when extended_promoters is
               active, so a promoter whose core lies in a neighbouring region
               but whose extended bounds reach into this one is still found.
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub _get_regfeats_by_coords {
  my ($self, $c, $s, $e) = @_;

  my $parser = $self->parser();

  my $source_chr = $self->get_source_chr_name($c, 'regfeat', $self->valid_chromosomes);

  # pad the query on both sides so a promoter whose core sits in a
  # neighbouring bin, but whose --extended_promoters reach crosses into this
  # one, is not invisible to this bin's seek. See $EXTENDED_PROMOTER_MARGIN.
  my $margin = $self->{extended_promoters} ? $EXTENDED_PROMOTER_MARGIN : 0;

  return [] unless $parser->seek($source_chr, $s - 1 - $margin, $e + 1 + $margin);

  my $include = $self->include_feature_types;
  my @features;

  $parser->next();

  while($parser->{record} && $parser->get_start <= $e + $margin) {
    my $type = $parser->get_type;

    if($include->{$type}) {
      my $rf = $self->_record_to_regfeat($c, $type);
      push @features, $rf if $rf;
    }
    else {
      $self->warning_msg(
        "WARNING: Ignoring '$type' feature type from ".$self->file.
        "; not a supported regulatory feature type\n"
      ) unless $self->{_warned_type}->{$type}++;
    }

    $parser->next();
  }

  return \@features;
}


=head2 _record_to_regfeat

  Arg 1      : string $chr
  Arg 2      : string $type
  Example    : $rf = $as->_record_to_regfeat($chr, $type);
  Description: Converts the parser's current record into a
               Bio::EnsEMBL::Funcgen::RegulatoryFeature carrying only the fields
               that survive into the VEP cache. The md5 of the raw record is
               attached so merge_features() can deduplicate features returned
               from more than one cache region.
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : none
  Caller     : _get_regfeats_by_coords()
  Status     : Stable

=cut

sub _record_to_regfeat {
  my ($self, $chr, $type) = @_;

  my $parser = $self->parser;
  my $attributes = $parser->get_attributes || {};

  my ($start, $end) = ($parser->get_start, $parser->get_end);

  # --extended_promoters widens promoters to their extended bounds; the
  # extension is always additive and only promoters carry these attributes
  # (motif GFFs have none, and is_motif sources never set extended_promoters)
  if($self->{extended_promoters} && $HAS_EXTENDED_BOUNDS{$type}) {
    my $xs = $attributes->{extended_start};
    my $xe = $attributes->{extended_end};
    $start = $xs if defined($xs) && $xs < $start;
    $end   = $xe if defined($xe) && $xe > $end;
  }

  my $slice = $self->_get_or_fake_slice($chr);

  my $strand = $parser->get_strand;
  $strand = 0 unless defined($strand);

  my $stable_id = $self->_record_get_id($attributes, $chr, $start, $end, $type);

  my $feature = $self->{is_motif}
    ? $self->_build_motif_feature($attributes, $slice, $start, $end, $strand, $stable_id)
    : $self->_build_regulatory_feature($slice, $start, $end, $strand, $stable_id, $type);

  # deduplication key: md5 of the raw record. Features straddling a cache
  # region boundary are read once per region and must collapse to one.
  $feature->{md5} = $self->_record_md5;

  $self->_check_id_collision($feature);

  return $feature;
}


=head2 _build_regulatory_feature

  Description: Builds a Bio::EnsEMBL::Funcgen::RegulatoryFeature from the current
               record.
  Returntype : Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Caller     : _record_to_regfeat()
  Status     : Stable

=cut

sub _build_regulatory_feature {
  my ($self, $slice, $start, $end, $strand, $stable_id, $type) = @_;

  my $rf = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new_fast({
    start             => $start,
    end               => $end,
    strand            => $strand,
    slice             => $slice,
    stable_id         => $stable_id,
    dbID              => ++$self->{_rf_dbID},
    feature_type      => $type,
    _vep_feature_type => 'RegulatoryFeature',
  });

  return $rf;
}


=head2 _build_motif_feature

  Description: Builds a Bio::EnsEMBL::Funcgen::MotifFeature from the current
               record. Carries a BindingMatrix stub with stable_id, length (the
               motif span) and the transcription factor complexes; elements is
               deliberately left unset, so OutputFactory omits HIGH_INF_POS and
               MOTIF_SCORE_CHANGE (which need the position weight matrix). A
               future matrices= source fills elements in and both fields light up
               with no further change.
  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature
  Caller     : _record_to_regfeat()
  Status     : Stable

=cut

sub _build_motif_feature {
  my ($self, $attributes, $slice, $start, $end, $strand, $stable_id) = @_;

  # transcription_factor is a comma-separated list; model each as a complex with
  # a display_name, which is all OutputFactory reads
  my @tfcs =
    map { { display_name => $_ } }
    split(/,/, $attributes->{transcription_factor} || '');

  my $matrix = bless {
    stable_id                                 => $attributes->{binding_matrix_id},
    length                                    => ($end - $start + 1),
    associated_transcription_factor_complexes => \@tfcs,
  }, 'Bio::EnsEMBL::Funcgen::BindingMatrix';

  return Bio::EnsEMBL::Funcgen::MotifFeature->new_fast({
    start             => $start,
    end               => $end,
    strand            => $strand,
    slice             => $slice,
    stable_id         => $stable_id,
    dbID              => ++$self->{_rf_dbID},
    binding_matrix    => $matrix,
    _vep_feature_type => 'MotifFeature',
  });
}


=head2 _check_id_collision

  Arg 1      : Bio::EnsEMBL::Funcgen::RegulatoryFeature $rf
  Example    : $as->_check_id_collision($rf);
  Description: Warns if two distinct records share a stable identifier.
               Bio::EnsEMBL::Variation::VariationFeature keys its
               regulatory_feature_variations hash on the feature's stable id, so
               a collision silently discards one of the annotations - last write
               wins. Distinguished from the legitimate case of one feature being
               read from two adjacent cache regions by comparing the record md5:
               the same record read twice has the same md5, two different records
               do not.
  Returntype : none
  Exceptions : none
  Caller     : _record_to_regfeat()
  Status     : Stable

=cut

sub _check_id_collision {
  my ($self, $rf) = @_;

  my $id  = $rf->{stable_id};
  my $md5 = $rf->{md5};
  return unless defined($id) && defined($md5);

  my $seen = $self->{_seen_ids} ||= {};

  if(exists($seen->{$id}) && $seen->{$id} ne $md5) {
    $self->warning_msg(
      "WARNING: Duplicate regulatory feature identifier '$id' in ".$self->file.
      "; identifiers must be unique or annotations will be lost"
    );
  }
  else {
    $seen->{$id} = $md5;
  }

  return;
}


=head2 _get_or_fake_slice

  Arg 1      : string $chr
  Example    : $slice = $as->_get_or_fake_slice($chr);
  Description: Gets a slice for the given chromosome, fabricating one if none is
               available. Regulatory consequences are purely coordinate-based -
               within_regulatory_feature reduces to a start/end overlap - so no
               sequence is required, but BaseVEP::get_slice returns undef when
               there is no sequence source at all. Without this fallback
               _record_to_regfeat would return undef and every feature would be
               dropped silently.

               Reached only when running --cache --offline against a cache with
               no FASTA alongside it: there is then no core adaptor and no
               fasta_db, so no slice can be built. The primary use case -
               --gff transcripts, which require --fasta, plus --regulatory_gff -
               always has a real slice and never reaches this. Removing it in
               favour of requiring --fasta would therefore only affect the cache
               case, at the cost of a genome download for sequence that is never
               read.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : _record_to_regfeat()
  Status     : Stable

=cut

sub _get_or_fake_slice {
  my ($self, $chr) = @_;

  return $self->{_slice_for}->{$chr} ||= do {
    my $slice = $self->get_slice($chr);

    unless($slice) {
      $slice = Bio::EnsEMBL::Slice->new_fast({
        coord_system      => $self->{_coord_system} ||= Bio::EnsEMBL::CoordSystem->new(
          -NAME => 'chromosome', -RANK => 1
        ),
        start             => 1,
        end               => $FAKE_SLICE_LENGTH,
        seq_region_name   => $chr,
        seq_region_length => $FAKE_SLICE_LENGTH,
      });
      $slice->{is_fake} = 1;
    }

    $slice;
  };
}


=head2 _record_get_id

  Arg 1      : hashref $attributes
  Arg 2      : string $chr
  Arg 3      : int $start
  Arg 4      : int $end
  Arg 5      : string $type
  Example    : $id = $as->_record_get_id($attributes, $chr, $start, $end, $type);
  Description: Gets a stable identifier for a record, falling back to a
               deterministic coordinate-derived value so that a third-party GFF
               without ID or Name attributes still yields named features.

               The fallback is derived rather than counted because the value
               reaches the Feature column of the output: a counter would depend
               on read order, and a feature straddling a cache region boundary is
               built once per region, so the two builds must agree.
  Returntype : string
  Exceptions : none
  Caller     : _record_to_regfeat()
  Status     : Stable

=cut

sub _record_get_id {
  my ($self, $attributes, $chr, $start, $end, $type) = @_;

  for my $key(qw(ID Name)) {
    my $id = $attributes->{$key};
    return $id if defined($id) && length($id);
  }

  return sprintf('%s:%d-%d:%s', $chr, $start, $end, $type);
}


=head2 _record_md5

  Example    : $md5 = $as->_record_md5();
  Description: md5 of the parser's current raw record, used as the deduplication
               key in merge_features(). Identical for the same record however
               many cache regions it is read from.
  Returntype : string
  Exceptions : none
  Caller     : _record_to_regfeat()
  Status     : Stable

=cut

sub _record_md5 {
  my $self = shift;
  my $parser = $self->parser;

  my $raw = $parser->{current_block};
  $raw = join("\t", @{$parser->{record} || []}) unless defined($raw);

  require Digest::MD5;
  return Digest::MD5::md5_hex($raw);
}


=head2 merge_features

  Arg 1      : arrayref of Bio::EnsEMBL::Funcgen::RegulatoryFeature $features
  Example    : $merged = $as->merge_features($features);
  Description: Deduplicates features on the md5 of their source record.
               Overrides AnnotationType::RegFeat::merge_features, which keys on
               dbID; dbIDs here are assigned by an incrementing counter, so the
               same feature read from two adjacent cache regions would be given
               two different dbIDs and survive as a duplicate.
  Returntype : arrayref of Bio::EnsEMBL::Funcgen::RegulatoryFeature
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

sub merge_features {
  my ($self, $features) = @_;

  my (@return, %seen);

  foreach my $f(@$features) {
    my $key = defined($f->{md5}) ? $f->{md5} : $f->{stable_id};
    next if $key && $seen{$key}++;
    push @return, $f;
  }

  return \@return;
}


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $as->annotate_InputBuffer($ib);
  Description: Delegates to the regulatory implementation. Required because
               AnnotationSource::File does not provide one that creates
               RegulatoryFeatureVariation objects.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub annotate_InputBuffer {
  my ($self, $buffer) = @_;
  return Bio::EnsEMBL::VEP::AnnotationType::RegFeat::annotate_InputBuffer($self, $buffer);
}


=head2 info

  Example    : $info = $as->info()
  Description: Gets the info hashref for this annotation source.
  Returntype : hashref
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::BaseRunner
  Status     : Stable

=cut

sub info {
  my $self = shift;
  return $self->{info} ||= { regulatory_gff => $self->file };
}

1;
