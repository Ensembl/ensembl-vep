=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Core;

use strict;
use warnings;

use Scalar::Util qw(weaken);
use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript;

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper);

sub run {
  my $self = shift;

  my $vep_params = $self->get_vep_params();

  $vep_params->{$_} = 1 for qw(
    protein
    xref_refseq
    domains
    uniprot
  );

  # core doesn't use BAM (yet?)
  delete $vep_params->{bam} if $vep_params->{bam};

  my $config = Bio::EnsEMBL::VEP::Config->new($vep_params);
  $config->{_registry} = 'Bio::EnsEMBL::Registry';
  my $region_size = $self->param('region_size');

  my $hive_dbc = $self->dbc;
  $hive_dbc->disconnect_if_idle() if defined $hive_dbc;

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $config,
    cache_region_size => $region_size,
  });

  $config->param($_, 0) for qw(sift polyphen);

  my $cache = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
    config => $config,
    cache_region_size => $region_size,
    dir => $self->get_cache_dir($vep_params)
  });

  $self->dump_chrs($as, $cache);

  $self->dump_info($as, $self->get_cache_dir($vep_params));
  
  Bio::EnsEMBL::Registry->get_DBAdaptor($vep_params->{species},'core')->dbc()->disconnect_if_idle();

  return;
}

sub dump_info {
  my ($self, $as, $dir, $extra) = @_;

  my $info_file = $dir.'/info.txt_core';
  return if -e $info_file;

  open OUT, ">$info_file" or die "Could not write to info file $info_file\n";

  print OUT "$_\t".$self->param($_)."\n" for qw(species assembly sift polyphen);

  print OUT "$_\t".$extra->{$_}."\n" for keys %{$extra || {}};

  my $info = $as->info;
  print OUT "source_$_\t".$info->{$_}."\n" for keys %$info;

  close OUT;
}

sub get_dumpable_object {
  my ($self, $as, $sr, $chr, $s) = @_;

  my $obj = {
    $chr => $as->merge_features([
      map {
        $as->apply_edits($_);
        $as->lazy_load_transcript($_);
        $self->clean_transcript($_);
        $_
      }
      @{$as->get_features_by_regions_uncached([[$sr, $s]], 1)}
    ])
  };

  for (my $i = 0; $i < scalar(@{$obj->{$chr}}); $i++) {
    if(my $tr = $obj->{$chr}->[$i]) {
      delete $tr->{slice}->{adaptor};
      delete $tr->{slice}->{coord_system}->{adaptor};

      # clean introns, they may have a rev-strand slice attached
      if(my ($intron) = @{$tr->{_variation_effect_feature_cache}->{introns}}) {
        delete $intron->{slice}->{adaptor};
        delete $intron->{slice}->{coord_system}->{adaptor};
      }
    }
  }

  return $obj;
}

sub post_dump {
  my ($self, $obj, $as, $chr) = @_;

  if(my $tr = $obj->{$chr}->[0]) {
    $tr->{slice}->{adaptor} = $as->get_adaptor('core', 'slice');
    $tr->{slice}->{coord_system}->{adaptor} = $as->get_adaptor('core', 'coordsytem');

    if(my ($intron) = @{$tr->{_variation_effect_feature_cache}->{introns}}) {
      $intron->{slice}->{adaptor} = $as->get_adaptor('core', 'slice');;
      $intron->{slice}->{coord_system}->{adaptor} = $as->get_adaptor('core', 'slice');;
    }
  }
}

sub clean_transcript {
  my $self = shift;
  my $tr = shift;

  # delete keys
  foreach my $key(qw(
    display_xref
    external_db
    external_display_name
    external_name
    external_status
    created_date
    status
    description
    edits_enabled
    modified_date
    dbentries
    is_current
    analysis
    transcript_mapper
    adaptor
  )) {
    delete $tr->{$key} if defined($tr->{$key});
  }

  # clean attributes
  if(defined($tr->{attributes})) {
    my @new_atts;
    my %keep = map {$_ => 1} qw(
      gencode_basic
      miRNA
      ncRNA
      cds_start_NF
      cds_end_NF
      TSL
      appris
      rseq_mrna_match
      rseq_mrna_nonmatch
      rseq_5p_mismatch
      rseq_cds_mismatch
      rseq_3p_mismatch
      rseq_nctran_mismatch
      rseq_no_comparison
      rseq_ens_match_wt
      rseq_ens_match_cds
      rseq_ens_no_match
      enst_refseq_compare
      _rna_edit
      MANE_Select
      MANE_Plus_Clinical 
   );

    foreach my $att(@{$tr->{attributes}}) {
      delete $att->{description};
      push @new_atts, $att if defined($keep{$att->{code}});
    }
    $tr->{attributes} = \@new_atts;
  }

  # clean exons
  foreach my $exon(@{$tr->{_trans_exon_array}}) {
    delete $exon->{$_} for qw(
      adaptor
      created_date
      modified_date
      is_current
      version
      is_constitutive
      _seq_cache
      dbID
      slice
    );
  }

  # clean the translation
  if($tr->translation) {

    # sometimes the translation points to a different transcript?
    $tr->{translation}->{transcript} = $tr;
    weaken($tr->{translation}->{transcript});

    for my $key(qw(attributes protein_features created_date modified_date dbentries adaptor)) {
      delete $tr->translation->{$key};
    }
  }

  # clean gene
  if(my $gene = $tr->{_gene}) {
    my %keep = map {$_ => 1} qw{start end strand stable_id};
    map {delete $gene->{$_}} grep {!$keep{$_}} keys %{$gene};
  }
}

1;
