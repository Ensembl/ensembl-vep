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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript - database transcript annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
  config => $config,
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Database-based annotation source for transcript data.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript;

use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Database
  Bio::EnsEMBL::VEP::AnnotationType::Transcript
);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config $config,
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript
  Exceptions : warns if --nearest set
  Caller     : AnnotationSourceAdaptor
  Status     : Stable

=cut

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
    merged
    sift
    polyphen
    everything
    xref_refseq
    protein
    uniprot
    domains
    use_transcript_ref
    no_prefetch
  )]);

  if($self->param('nearest')) {
    $self->warning_msg("WARNING: --nearest is not currently compatible with database annotation type");
  }

  $self->{source_type} = ($self->{core_type} || '') eq 'otherfeatures' ? 'refseq' : 'ensembl';

  $self->check_sift_polyphen();

  return $self;
}


=head2 check_sift_polyphen

  Example    : $ok = $as->check_sift_polyphen();
  Description: Gets user-set SIFT/PolyPhen parameters and checks vs
               availability in the database. If using "safe" mode (REST, web)
               or --everything, params are disabled if unavailable. Otherwise,
               this method will throw.
  Returntype : bool
  Exceptions : throws if configured tool not available
  Caller     : new()
  Status     : Stable

=cut

sub check_sift_polyphen {
  my $self = shift;

  my $info = $self->info;

  foreach my $tool(grep {$self->{lc($_)}} qw(SIFT PolyPhen)) {
    my $lc_tool = lc($tool);

    unless($info->{$lc_tool}) {

      # dont die if user set "everything" param on a species with no SIFT/PolyPhen
      if($self->{everything} || $self->param('safe')) {
        $self->status_msg("INFO: disabling $tool");
        $self->param($lc_tool, 0);
        $self->{$lc_tool} = 0;
      }

      else {
        throw("ERROR: $tool not available\n");
      }
    }
  }

  return 1;
}


=head2 assembly

  Example    : $assembly = $as->assembly();
  Description: Get the default assembly version for the underlying database.
  Returntype : string
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub assembly {
  my $self = shift;
  return $self->{assembly} ||= $self->get_database_assembly;
}


=head2 get_features_by_regions_uncached

  Arg 1      : arrayref $regions
  Example    : $trs = $as->get_features_by_regions_uncached($regions)
  Description: Gets all transcripts overlapping the given set of regions. See
               Bio::EnsEMBL::VEP::AnnotationSource::get_all_regions_by_InputBuffer()
               for information about regions.
  Returntype : arrayref of Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer(), DumpVEP pipeline
  Status     : Stable

=cut

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;
  my $chr_is_seq_region = shift;

  my $cache = $self->cache;
  my @return;

  my $pfa = $self->get_adaptor('variation', 'PhenotypeFeature');

  my $cache_region_size = $self->{cache_region_size};

  foreach my $region(@{$regions}) {
    my ($c, $region_start) = @$region;

    my $slice = $self->get_slice($c, $chr_is_seq_region);

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
        next if $tr->analysis() && $tr->analysis()->logic_name() eq 'estgene';
        # in human and mouse otherfeatures DB, there may be duplicate genes
        # skip those from analysis refseq_human_import and refseq_mouse_import
        next if $self->{core_type} eq 'otherfeatures' && $self->assembly !~ /GRCh37/i && $tr->analysis && $tr->analysis->logic_name =~ /^refseq_[a-z]+_import$/;

	## Due to the inclusion of the new RefSeq transcript set (mapped from 38) into the 37 otherfeatures database,
	## older, lower quality transcripts have been removed from the cache files. To do this, we filter out all transcripts
	## with the analysis logic_name 'refseq_import' or 'refseq_human_import'. 
	## All new transcripts have the logic name 'refseq_import_grch38'
	next if $self->{core_type} eq 'otherfeatures' && $self->assembly =~ /GRCh37/i && $tr->analysis && $tr->analysis->logic_name =~ /^refseq.+import$/;
	if($self->{core_type} eq 'otherfeatures' && defined($tr->display_xref)){
	         $tr->{stable_id} = $tr->display_xref->{display_id};
        }
        
        $tr->{_gene_stable_id} = $gene_stable_id;
        $tr->{_gene} = $gene;
        $self->prefetch_gene_ids($tr);

        # flag transcript source if using "merged" (ie both core and otherfeatures DBs)
        $tr->{_source_cache} = $self->{core_type} eq 'otherfeatures' ? 'RefSeq' : 'Ensembl' if $self->{merged};

        # indicate if canonical
        $tr->{is_canonical} = 1 if defined $canonical_tr_id and $tr->dbID eq $canonical_tr_id;

        # indicate phenotype
        $tr->{_gene_phenotype} = $gene_has_phenotype;

        push @features, $tr;
      }
    }

    $cache->{$c}->{$region_start} = \@features;

    push @return, @features;
  }

  return \@return;
}


=head2 lazy_load_transcript

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Example    : $tr = $as->lazy_load_transcript($tr)
  Description: Prefetches all transcript data (sequence, identifiers etc).
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub lazy_load_transcript {
  my ($self, $tr) = @_;
  
  unless($tr->{_vep_lazy_loaded}) {
    $self->prefetch_transcript_data($tr) unless $self->{no_prefetch};
    $tr->{_vep_lazy_loaded} = 1;
  }

  return $tr;
}


=head2 prefetch_transcript_data

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Example    : $tr = $as->prefetch_transcript_data($tr)
  Description: Prefetches all transcript data; data are cached on
               key named "_variation_effect_feature_cache" on transcript
               object. Includes:
                - introns
                - exons (sorted by chromosome position)
                - translateable sequence
                - transcript mapper
                - 3' UTR sequence
                - codon table
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : lazy_load_transcript()
  Status     : Stable

=cut

sub prefetch_transcript_data {
  my $self = shift;
  my $tr = shift;

  my $vep_cache = $tr->{_variation_effect_feature_cache} ||= {};

  $vep_cache->{introns} ||= $tr->get_all_Introns;
  $vep_cache->{sorted_exons} ||= [sort {$a->start <=> $b->start} @{$tr->get_all_Exons}];
  $vep_cache->{translateable_seq} ||= $tr->translateable_seq;
  $vep_cache->{mapper} ||= $tr->get_TranscriptMapper;

  # UTRs
  my $transferred = $tr->transfer($tr->feature_Slice());
  $vep_cache->{$_.'_prime_utr'} = $self->_fetch_utr($_, $transferred) for qw(five three);

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


=head2 _fetch_utr

  Arg 1      : string $five_three ('five' or 'three')
  Arg 2      : Bio::EnsEMBL::Transcript $tr
  Example    : $utr = $as->_fetch_utr('five', $tr)
  Description: Get Bio::Seq object representing either the five- or
               three-prime UTR. Uses a cache so that multiple transcripts
               with the same UTR sequence will share a reference to the
               same Bio::Seq object.
  Returntype : Bio::Seq
  Exceptions : warns if fail eval to fetch UTR
  Caller     : prefetch_transcript_data()
  Status     : Stable

=cut

sub _fetch_utr {
  my $self = shift;
  my $five_three = shift;
  my $tr = shift;

  my $key = $five_three.'_prime_utr';
  my $utr;

  eval {
    if(my $this_utr = $tr->$key()) {
      my $md5 = md5_hex($this_utr->seq);

      my $utr_obj;
      my $cache = $self->{'_utr_cache'} ||= [];

      unless(($utr_obj) = map {$_->{obj}} grep {$_->{md5} eq $md5} @$cache) {
        $utr_obj = $this_utr;
        push @$cache, {md5 => $md5, obj => $utr_obj};
        shift @$cache if scalar @$cache > 50;
      }

      $utr = $utr_obj;
    }
  };
  if($@) {
    warn "Problem getting $five_three prime UTR:".$@;
  }

  return $utr;
}


=head2 prefetch_translation_data

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Arg 2      : (optional) Bio::EnsEMBL::Translation $tl
  Example    : $tr = $as->prefetch_translation_data($tr, $tl)
  Description: Prefetches all translation (protein) data; data are cached on
               key named "_variation_effect_feature_cache" on transcript
               object. Includes:
                - peptide sequence
                - protein domains
                - SeqEdits (e.g. selenocysteines)
                - SIFT and PolyPhen prediction matrices
                - identifiers
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : prefetch_transcript_data()
  Status     : Stable

=cut

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
  $vep_cache->{protein_function_predictions} = $self->get_sift_polyphen($vep_cache) if($self->{sift} || $self->{polyphen});

  $self->prefetch_translation_ids($tr, $tl);

  return $tr;
}


=head2 get_sift_polyphen

  Arg 1      : hashref $vep_cache
  Example    : $data = $as->get_sift_polyphen($tr->{_variation_effect_feature_cache})
  Description: Gets SIFT/PolyPhen predictions. Uses an internal cache so that transcripts
               with identical translations can share a reference to the same prediction
               matrix.
  Returntype : hashref
  Exceptions : none
  Caller     : prefetch_translation_data()
  Status     : Stable

=cut

sub get_sift_polyphen {
  my $self = shift;
  my $vep_cache = shift;

  my $md5 = md5_hex($vep_cache->{peptide});
  my $cache = $self->{_sift_polyphen_cache} ||= [];
  my $data;

  unless(($data) = map {$_->{data}} grep {$_->{md5} eq $md5} @$cache) {
    if(my $pfpma = $self->get_adaptor('variation', 'ProteinFunctionPredictionMatrix')) {
      foreach my $a(qw(sift polyphen_humvar polyphen_humdiv)) {
        next unless defined($self->{(split "_", $a)[0]});
        $data->{$a} ||= $pfpma->fetch_by_analysis_translation_md5($a, $md5);
        delete $data->{$a}->{adaptor};
      }
    }

    push @$cache, {md5 => $md5, data => $data};
    shift @$cache if scalar @$cache > 50;
  }

  return $data;
}


=head2 prefetch_gene_ids

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Example    : $tr = $as->prefetch_gene_ids($tr)
  Description: Prefetches gene identifiers for this transcript:
               - gene symbol
               - gene symbol source
               - HGNC ID
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub prefetch_gene_ids {
  my $self = shift;
  my $tr = shift;

  $tr->{_gene} ||= $tr->get_Gene();
  if( $self->{core_type} eq 'otherfeatures'){
    #Pulls correct gene symbol for RefSeq using xrefs
    my @entries = grep {$_->{dbname} eq 'EntrezGene'} @{$tr->get_Gene()->get_all_DBEntries};
    if(scalar @entries eq 1)
    {
      my $xref_obj = $tr->{_gene}->display_xref;
      $tr->get_Gene()->stable_id($entries[0]->{primary_id});
      $tr->{_gene_symbol}  = $xref_obj ? $xref_obj->display_id : $entries[0]->{display_id};
      $tr->{_gene_symbol_source} = $entries[0]->{dbname};
      $tr->{_gene_symbol_id} = $entries[0]->{primary_id};
      $tr->{_gene_hgnc_id} = $entries[0]->{primary_id} if $entries[0]->{dbname} eq 'HGNC';
      $tr->{_gene_stable_id} = $entries[0]->{primary_id};
    }
  }
  else {
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

  }

  return $tr;
}


=head2 prefetch_transcript_ids

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Example    : $tr = $as->prefetch_transcript_ids($tr)
  Description: Prefetches transcript identifiers for this transcript:
               - CCDS
               - RefSeq (via xref)
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : prefetch_transcript_data()
  Status     : Stable

=cut

sub prefetch_transcript_ids {
  my $self = shift;
  my $tr = shift;

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


=head2 prefetch_translation_ids

  Arg 1      : Bio::EnsEMBL::Transcript $tr
  Arg 2      : (optional) Bio::EnsEMBL::Translation $tl
  Example    : $tr = $as->prefetch_translation_ids($tr)
  Description: Prefetches translation identifiers for this transcript:
               - Uniprot (SWISSPROT, TrEMBL, Uniparc)
               - ENSP
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : prefetch_translation_data()
  Status     : Stable

=cut

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
    
    $tr->{_uniprot_isoform} = '-';
    @entries = grep {$_->database eq 'Uniprot_isoform'} @{$tl->get_all_DBEntries};
    if(scalar @entries) {
      $tr->{_uniprot_isoform} = join ",", map {$_->display_id} @entries;
    }

  }

  # Ensembl protein ID
  if($self->{protein}) {
    $tr->{_protein} = $tl->stable_id;
    # With the new RefSeq otherfeatures database, RefSeq identifiers can be accessed through get_all_DB_Entries.
    # Identifiers are accessed this way and relavent objects updated.
    my @entries = grep {$_->{dbname} eq 'GenBank'} @{$tl->get_all_DBEntries};
    if(scalar @entries == 1)
    {
      $tr->{_protein} = $entries[0]->{primary_id};
      $tl->{stable_id} = $entries[0]->{primary_id};
    }
  }

  return $tr;
}


=head2 info

  Example    : $info = $as->info()
  Description: Gets the info hashref for this annotation source. Contains
               version information for:
                - assembly name
                - GENCODE version
                - genebuild date
                - SIFT and PolyPhen versions
                - RefSeq GFF version (where applicable)
  Returntype : hashref
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::BaseRunner
  Status     : Stable

=cut

sub info {
  my $self = shift;

  if(!exists($self->{info})) {
    my %info;

    # core source versions
    if(my $core_mca = $self->get_adaptor('core', 'metacontainer')) {
      foreach my $meta_key(qw(assembly.name gencode.version genebuild.initial_release_date genebuild.version)) {
        my $version = $core_mca->list_value_by_key($meta_key);

        my $new_key = $meta_key;
        $new_key =~ s/\..+//;
        $info{$new_key} = $version->[0] if defined($version) && scalar @$version;
      }
    }

    # sift/polyphen versions
    if(my $var_mca = $self->get_adaptor('variation', 'metacontainer')) {
      foreach my $tool(grep {$self->{$_}} qw(sift polyphen)) {
        my $sth = $var_mca->db->dbc->prepare(qq{
          SELECT meta_value
          FROM meta
          WHERE meta_key = ?
        });
        $sth->execute($tool.'_version');

        my $version;
        $sth->bind_columns(\$version);
        $sth->fetch();
        $sth->finish();
        
        $info{$tool} = $version if defined($version);
      }
    }

    # refseq
    if($self->{source_type} eq 'refseq') {

      my $logic_name = 'refseq_import';
      if($info{'assembly'} =~ /GRCh37/) {
        $logic_name = 'refseq_import_grch38';
      }

      if(my $refseq_mca = $self->get_adaptor('otherfeatures', 'metacontainer')) {
        my $sth = $refseq_mca->db->dbc->prepare(qq{
          SELECT CONCAT(db_version, IF(db_file IS NOT NULL, concat(' - ', db_file), ''))
          FROM analysis
          WHERE logic_name = '$logic_name'
        });
        $sth->execute;

        my $version;
        $sth->bind_columns(\$version);
        $sth->fetch;
        $sth->finish;

        $info{refseq} = $version if defined($version);
      }
    }

    $self->{info} = \%info;
  }

  return $self->{info};
}

1;
