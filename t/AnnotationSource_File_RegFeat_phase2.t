# Copyright [2016-2026] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Test::Warnings qw(warning :no_end_test);
use FindBin qw($Bin);

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

my $reg      = $test_cfg->{regulatory_gff};
my $motifs   = $test_cfg->{regulatory_gff_motifs};
my $activity = $test_cfg->{regulatory_gff_activity};

use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat');
use_ok('Bio::EnsEMBL::VEP::AnnotationSourceAdaptor');
use_ok('Bio::EnsEMBL::VEP::Config');

sub cfg { Bio::EnsEMBL::VEP::Config->new({ %{$test_cfg->base_testing_cfg}, @_ }) }


## MOTIF FEATURE SOURCE
#######################

{
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => cfg(regulatory => 1), file => $motifs, motif => 1
  });
  is(ref($as), 'Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat', 'motif source class');
  is($as->short_name, 'MotifFeatures', 'motif short_name');

  my $feats = $as->_get_regfeats_by_coords(21, 25585700, 25585760);
  ok($feats && @$feats, 'motif features read');

  # class gate: MotifFeature consequences depend on the Perl class via isa()
  is(ref($feats->[0]), 'Bio::EnsEMBL::Funcgen::MotifFeature', 'built as MotifFeature');
  is($feats->[0]->{_vep_feature_type}, 'MotifFeature', '_vep_feature_type is MotifFeature');

  my ($mf) = grep {$_->stable_id eq 'ENSM21_TFBS1'} @$feats;
  ok($mf, 'motif by stable_id');
  is($mf->strand, -1, 'motif strand');

  # BindingMatrix stub: identifier, TF complexes, length (from feature span);
  # no elements, so the PWM-dependent fields are unavailable downstream
  my $bm = $mf->get_BindingMatrix;
  ok($bm, 'has BindingMatrix');
  is($bm->stable_id, 'ENSPFM0001', 'matrix stable_id from binding_matrix_id');
  is($bm->{length}, 21, 'matrix length from feature span');
  is_deeply(
    [ map {$_->{display_name}} @{$bm->{associated_transcription_factor_complexes}} ],
    ['CTCF', 'FOO'],
    'transcription factors from attribute'
  );
  ok(!defined $bm->{elements}, 'BindingMatrix has no PWM elements');
}


## OUTPUT: TF_binding_site_variant, with PWM fields omitted
##########################################################

{
  use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
  use_ok('Bio::EnsEMBL::VEP::InputBuffer');
  use_ok('Bio::EnsEMBL::VEP::OutputFactory::Tab');

  my $cfg = cfg(regulatory => 1, fasta => $test_cfg->{fasta});
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => $cfg, file => $motifs, motif => 1
  });

  my $vcf = $test_cfg->create_input_file([['21', 25585730, 'in_motif', 'C', 'T']]);
  my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({config => $cfg, file => $vcf, valid_chromosomes => [21]});
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();
  $as->annotate_InputBuffer($ib);

  my ($vf) = @{$ib->buffer};
  my @mfvs = values %{$vf->{motif_feature_variations} || {}};
  is(scalar @mfvs, 1, 'one MotifFeatureVariation created');

  $_->{transcript_variations} ||= {} for @{$ib->buffer};

  my $tab = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({config => $cfg});
  my ($line) = grep {/MotifFeature/} @{$tab->get_all_lines_by_InputBuffer($ib)};
  ok($line, 'motif output line produced');
  like($line, qr/TF_binding_site_variant/, 'consequence is TF_binding_site_variant');
  like($line, qr/ENSPFM0001/, 'MOTIF_NAME present');
  like($line, qr/CTCF/, 'TRANSCRIPTION_FACTORS present');

  # HIGH_INF_POS and MOTIF_SCORE_CHANGE need the PWM, which the GFF lacks - they
  # must be omitted, not printed as a (misleading) value
  my $hash = $tab->get_all_output_hashes_by_InputBuffer($ib);
  my ($mh) = grep {($_->{Feature_type}||'') eq 'MotifFeature'} @$hash;
  ok($mh, 'motif output hash');
  ok(!exists $mh->{HIGH_INF_POS}, 'HIGH_INF_POS omitted without PWM');
  ok(!exists $mh->{MOTIF_SCORE_CHANGE}, 'MOTIF_SCORE_CHANGE omitted without PWM');
}


## --regulatory_gff key=value PARSING (via AnnotationSourceAdaptor)
###################################################################

# bare filename still works
{
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({
    config => cfg(regulatory_gff => $reg)
  });
  my $sources = $asa->get_all_regulatory_gff;
  is(scalar @$sources, 1, 'bare filename -> one source');
  is(ref($sources->[0]), 'Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat', 'bare -> RegFeat');
}

# file= plus motifs= plus emars= -> three sources of the right classes
{
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({
    config => cfg(regulatory_gff => "file=$reg,motifs=$motifs,emars=$reg")
  });
  my $sources = $asa->get_all_regulatory_gff;
  is(scalar @$sources, 3, 'file+motifs+emars -> three sources');
  is(ref($sources->[0]), 'Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat',  'main is RegFeat');
  is(ref($sources->[1]), 'Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat',  'emars is RegFeat');
  is($sources->[1]->short_name, 'EMARs', 'emars short_name');
  is(ref($sources->[2]), 'Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat', 'motifs source is File::RegFeat');
  is($sources->[2]->short_name, 'MotifFeatures', 'motifs source is in motif mode');
}

# unknown key rejected
{
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({
    config => cfg(regulatory_gff => "file=$reg,bogus=$reg")
  });
  throws_ok { $asa->get_all_regulatory_gff } qr/Unsupported --regulatory_gff key/, 'unknown key throws';
}

# file= required in key=value form
{
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({
    config => cfg(regulatory_gff => "motifs=$motifs")
  });
  throws_ok { $asa->get_all_regulatory_gff } qr/requires file=/, 'missing file= throws';
}

# matrices= reserved but not yet implemented
{
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({
    config => cfg(regulatory_gff => "file=$reg,matrices=$reg")
  });
  throws_ok { $asa->get_all_regulatory_gff } qr/matrices= is not yet supported/, 'matrices= rejected';
}


## EMAR type
############

{
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => cfg(regulatory => 1), file => $reg, short_name => 'EMARs'
  });
  my $feats = $as->_get_regfeats_by_coords(21, 28000000, 28002000);
  my ($emar) = grep {$_->stable_id eq 'ENSR21_EMAR1'} @$feats;
  ok($emar, 'EMAR feature read');
  is($emar->{feature_type}, 'EMAR', 'EMAR type carried to BIOTYPE');
  is($emar->{_vep_feature_type}, 'RegulatoryFeature', 'EMAR is a RegulatoryFeature');
}


## activity= join and real --cell_type support
##############################################

{
  # available cell types come from the activity header
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => cfg(regulatory => 1), file => $reg, activity => $activity
  });
  is_deeply($as->get_available_cell_types, ['GM12878', 'K562', 'HUVEC'], 'available cell types from activity header');
}

{
  # requested cell type is joined onto the matching feature
  my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => cfg(regulatory => 1, cell_type => ['GM12878']),
    file => $reg, activity => $activity
  });
  my $feats = $as->_get_regfeats_by_coords(21, 25585500, 25585900);
  my ($prom) = grep {$_->stable_id eq 'ENSR21_PROM1'} @$feats;
  ok($prom, 'feature with activity read');
  is($prom->{cell_types}->{GM12878}, 'ACTIVE', 'activity joined by stable_id');

  # a feature not in the activity table gets no cell types (the CTCF gap)
  my $ctcf = $as->_get_regfeats_by_coords(21, 25587690, 25587720);
  my ($c) = grep {$_->{feature_type} eq 'CTCF_binding_site'} @$ctcf;
  ok($c, 'CTCF feature read');
  ok(!%{$c->{cell_types} || {}}, 'feature absent from activity table has no cell types');
}

{
  # unknown cell type is rejected against the activity header
  throws_ok {
    Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
      config => cfg(regulatory => 1, cell_type => ['NOPE']),
      file => $reg, activity => $activity
    })
  } qr/unavailable/, 'unknown cell type rejected';
}

{
  # --cell_type without activity= still warns and is ignored (phase-1 behaviour)
  my $as;
  my $warned = warning {
    $as = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
      config => cfg(regulatory => 1, cell_type => ['GM12878']), file => $reg
    })
  };
  like("$warned", qr/--cell_type is ignored/, '--cell_type without activity warns');
  ok($as, 'source still built');
}

done_testing();
