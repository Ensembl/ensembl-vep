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

my $gff = $test_cfg->{regulatory_gff};

## BASIC TESTS
##############

use_ok('Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat');
use_ok('Bio::EnsEMBL::VEP::Config');

# base config: no cache, no database - a regulatory GFF needs neither
sub base_cfg {
  my %extra = @_;
  return Bio::EnsEMBL::VEP::Config->new({
    %{$test_cfg->base_testing_cfg},
    regulatory => 1,
    %extra,
  });
}

sub new_as {
  my %extra = @_;
  my $cfg = base_cfg(%extra);
  return Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => $cfg,
    file   => $gff,
  });
}

my $as = new_as();
ok($as, 'new is defined');
is(ref($as), 'Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat', 'check class');
is($as->file, $gff, 'file');
is($as->short_name, 'RegulatoryFeatures', 'short_name default');
is($as->{cache_region_size}, 1e6, 'cache_region_size');

throws_ok {
  Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({config => base_cfg()})
} qr/No file given/, 'no file throws';

# A GFF carries no epigenome activity, so --cell_type cannot be honoured. It is
# warned about and ignored rather than being an error, so the rest of the run
# still works; CELL_TYPE is simply empty for GFF-derived features.
{
  my $as_ct;
  my $warned = warning { $as_ct = new_as(cell_type => ['HUVEC']) };
  like("$warned", qr/--cell_type is ignored/, '--cell_type warns');
  ok($as_ct, '--cell_type does not prevent the source being created');
  ok(
    scalar @{$as_ct->_get_regfeats_by_coords(21, 25585000, 25586000)},
    '--cell_type does not prevent annotation'
  );
}


## PARSING AND OBJECT CONSTRUCTION
##################################

my $feats = $as->_get_regfeats_by_coords(21, 25585000, 25595000);
ok($feats && @$feats, '_get_regfeats_by_coords - got features');

# The class gate: consequence assignment blesses a dummy into the
# OverlapConsequence's feature_class and calls isa() on it. If these objects are
# not Funcgen::RegulatoryFeature there are no regulatory consequences AND no
# error, so this assertion guards a silent failure.
is(
  ref($feats->[0]),
  'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
  'features are Funcgen::RegulatoryFeature (guards silent consequence loss)'
);
ok(
  $feats->[0]->isa('Bio::EnsEMBL::Funcgen::RegulatoryFeature'),
  'feature isa Funcgen::RegulatoryFeature'
);
is($feats->[0]->{_vep_feature_type}, 'RegulatoryFeature', '_vep_feature_type set');

my %by_id = map {$_->stable_id => $_} @$feats;

# stable_id from the ID attribute
ok($by_id{ENSR21_PROM1}, 'stable_id from ID attribute');
is($by_id{ENSR21_PROM1}->{feature_type}, 'promoter', 'feature_type is GFF column 3');
is($by_id{ENSR21_PROM1}->{start}, 25585500, 'core start by default');
is($by_id{ENSR21_PROM1}->{end},   25585900, 'core end by default');

# each supported type is read and carried through verbatim
is($by_id{ENSR21_CTCF1}->{feature_type}, 'CTCF_binding_site', 'CTCF_binding_site type');
is($by_id{ENSR21_ENH1}->{feature_type},  'enhancer',          'enhancer type');
is($by_id{ENSR21_OCR1}->{feature_type},  'open_chromatin_region', 'open_chromatin_region type');

# strand is carried where present, defaulted to 0 where absent
is($by_id{ENSR21_CTCF1}->{strand}, 1, 'strand read from GFF');
is($by_id{ENSR21_ENH1}->{strand},  0, 'strand defaults to 0 when "."');

# identifier fallback chain: ID, then Name, then a deterministic derivation.
# Derived rather than counted: the value reaches the Feature column, and a
# straddling feature is built once per cache region, so both builds must agree.
ok($by_id{ENSR21_NAMEONLY}, 'stable_id falls back to Name attribute');
ok(
  (grep {$_->stable_id =~ /^21:25593000-25593100:enhancer$/} @$feats),
  'stable_id falls back to a deterministic coordinate-derived value'
);

# unsupported column 3 values are skipped, not fatal
ok(
  !(grep {$_->{feature_type} eq 'pseudogene'} @$feats),
  'unsupported feature type is skipped'
);


## EXTENDED PROMOTER BOUNDS
###########################

my $as_ext = new_as(extended_promoters => 1);
my $ext = $as_ext->_get_regfeats_by_coords(21, 25585000, 25595000);
my %ext_by_id = map {$_->stable_id => $_} @$ext;

is($ext_by_id{ENSR21_PROM1}->{start}, 25585000, '--extended_promoters widens start');
is($ext_by_id{ENSR21_PROM1}->{end},   25586200, '--extended_promoters widens end');

# the flag is promoter-only; nothing else carries extended attributes
is($ext_by_id{ENSR21_ENH1}->{start}, 25587740, '--extended_promoters leaves enhancer start');
is($ext_by_id{ENSR21_ENH1}->{end},   25587800, '--extended_promoters leaves enhancer end');


## DEDUPLICATION ACROSS CACHE REGION BOUNDARIES
###############################################

# ENSR21_STRADDLE spans the 26Mb boundary, so it is read once per adjacent cache
# region. dbIDs come from an incrementing counter and would differ between the
# two reads, so dedup must key on the record md5 instead.
{
  my $as_s = new_as();
  my $a = $as_s->_get_regfeats_by_coords(21, 25999000, 26000000);
  my $b = $as_s->_get_regfeats_by_coords(21, 26000001, 26001000);

  my @straddle = grep {$_->stable_id eq 'ENSR21_STRADDLE'} (@$a, @$b);
  is(scalar @straddle, 2, 'straddling feature is read once per region');
  isnt($straddle[0]->{dbID}, $straddle[1]->{dbID}, 'dbIDs differ between reads');
  is($straddle[0]->{md5}, $straddle[1]->{md5}, 'md5 is identical between reads');

  my $merged = $as_s->merge_features(\@straddle);
  is(scalar @$merged, 1, 'merge_features collapses the duplicate');

  # the dbID-keyed implementation this overrides would not
  my $merged2 = $as_s->merge_features([@$a, @$b]);
  my %seen;
  $seen{$_->stable_id}++ for @$merged2;
  ok(!(grep {$_ > 1} values %seen), 'no duplicate stable_ids survive merge_features');
}


## STABLE IDENTIFIER COLLISIONS
################################

# VariationFeature keys its regulatory_feature_variations hash on the feature's
# stable id, so two distinct features sharing one id silently discards an
# annotation - last write wins. Warn rather than fail, but do warn.
{
  my $as_dup = new_as();
  my $warned = warning { $as_dup->_get_regfeats_by_coords(21, 26999000, 27001000) };
  like(
    "$warned", qr/Duplicate regulatory feature identifier 'ENSR21_DUPID'/,
    'duplicate stable identifier warns'
  );
}

# A feature straddling a cache region boundary is legitimately built once per
# region with the same id. That must NOT be reported as a collision - the check
# distinguishes the two by comparing the record md5.
{
  my $as_str = new_as();
  my @warnings;
  {
    local $SIG{__WARN__} = sub { push @warnings, $_[0] };
    $as_str->_get_regfeats_by_coords(21, 25999000, 26000000);
    $as_str->_get_regfeats_by_coords(21, 26000001, 26001000);
  }
  ok(
    !(grep {/Duplicate regulatory feature identifier/} @warnings),
    'straddling feature re-read is not reported as a collision'
  );
}


## ANNOTATION AND CONSEQUENCES
##############################

use_ok('Bio::EnsEMBL::VEP::Parser::VCF');
use_ok('Bio::EnsEMBL::VEP::InputBuffer');

sub annotate {
  my ($source, $cfg, $file) = @_;
  my $p = Bio::EnsEMBL::VEP::Parser::VCF->new({
    config => $cfg, file => $file, valid_chromosomes => [21]
  });
  my $ib = Bio::EnsEMBL::VEP::InputBuffer->new({config => $cfg, parser => $p});
  $ib->next();
  $source->annotate_InputBuffer($ib);
  return $ib;
}

{
  my $cfg = base_cfg();
  my $source = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => $cfg, file => $gff
  });
  my $ib = annotate($source, $cfg, $test_cfg->{test_vcf});

  my ($vf) = grep {$_->{start} == 25585733} @{$ib->buffer};
  ok($vf, 'found variant overlapping the test promoter');

  my @rfvs = values %{$vf->{regulatory_feature_variations} || {}};
  is(scalar @rfvs, 1, 'one RegulatoryFeatureVariation created');
  is($rfvs[0]->feature->stable_id, 'ENSR21_PROM1', 'attached to the right feature');

  my %terms;
  for my $a (@{$rfvs[0]->get_all_alternate_RegulatoryFeatureVariationAlleles}) {
    $terms{$_->SO_term} = 1 for @{$a->get_all_OverlapConsequences};
  }
  ok($terms{regulatory_region_variant}, 'SNV gives regulatory_region_variant');

  # a variant in a CTCF site gets the same consequence - all regulatory types
  # collapse to regulatory_region_variant, BIOTYPE is the only discriminator
  my ($ctcf_vf) = grep {$_->{start} == 25587701} @{$ib->buffer};
  if($ctcf_vf) {
    my @c = values %{$ctcf_vf->{regulatory_feature_variations} || {}};
    my ($ctcf) = grep {$_->feature->stable_id eq 'ENSR21_CTCF1'} @c;
    ok($ctcf, 'CTCF binding site annotated');
    is($ctcf->feature->{feature_type}, 'CTCF_binding_site', 'BIOTYPE discriminates the type');
  }
}

# a variant inside the extended promoter but outside its core bounds is
# annotated only when --extended_promoters is set
{
  my $vcf = $test_cfg->create_input_file([
    ['21', 25585100, 'ext_only', 'A', 'G'],
  ]);

  for my $ext (0, 1) {
    my $cfg = base_cfg(extended_promoters => $ext);
    my $source = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
      config => $cfg, file => $gff
    });
    my $ib = annotate($source, $cfg, $vcf);
    my ($vf) = @{$ib->buffer};
    my $n = scalar values %{$vf->{regulatory_feature_variations} || {}};
    is($n, $ext ? 1 : 0, "extended-only variant annotated=".($ext?'yes':'no')." with extended_promoters=$ext");
  }
}

# structural variants pick up the ablation/amplification consequences for free,
# via AnnotationType::RegFeat's StructuralVariationOverlap branch
{
  my $cfg = base_cfg();
  my $source = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => $cfg, file => $gff
  });
  my $vcf = $test_cfg->create_input_file([
    ['21', 25585400, 'del', 'A', '<DEL>', '.', '.', 'SVTYPE=DEL;END=25586000'],
  ]);
  my $ib = annotate($source, $cfg, $vcf);
  my ($vf) = @{$ib->buffer};
  ok($vf, 'structural variant parsed');
  ok(
    exists($vf->{regulation_structural_variations}) ||
      scalar(values %{$vf->{regulatory_feature_variations} || {}}),
    'structural variant reaches the regulatory machinery'
  );
}


## CHROMOSOME SYNONYMS
######################

# The GFF uses Ensembl-style sequence names (1..22, X, Y) while input VCFs
# commonly use "chr1" style. Resolution happens in get_source_chr_name, called
# from _get_regfeats_by_coords with the tabix index's own seqnames.
{
  my $as_syn = new_as();

  ok($as_syn->chromosome_synonyms($test_cfg->{chr_synonyms}), 'load synonyms');
  is($as_syn->get_source_chr_name(21), 21, 'get_source_chr_name - exists, same');
  is($as_syn->get_source_chr_name('chr21'), 21, 'get_source_chr_name - strip chr');
  is($as_syn->get_source_chr_name('NC_000021.9'), 21, 'get_source_chr_name - synonym');

  # and end to end: a chr-prefixed input must still find the features
  my $feats = $as_syn->_get_regfeats_by_coords('chr21', 25585000, 25586000);
  ok(
    (grep {$_->stable_id eq 'ENSR21_PROM1'} @$feats),
    'chr-prefixed sequence name resolves to the GFF feature'
  );
}


## OUTPUT FORMATS
#################

# BIOTYPE carries the regulatory feature type and is the only field
# distinguishing an enhancer from a promoter, since every regulatory type
# collapses to regulatory_region_variant. Check it survives into both formats.
{
  use_ok('Bio::EnsEMBL::VEP::OutputFactory::Tab');
  use_ok('Bio::EnsEMBL::VEP::OutputFactory::JSON');

  # Output formatting resolves reference alleles, so it needs slices backed by
  # real sequence - a fabricated one is enough to compute consequences but not to
  # render them. Use the test FASTA, as a real run would use the cache or --fasta.
  my $cfg = base_cfg(fasta => $test_cfg->{fasta});
  my $source = Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat->new({
    config => $cfg, file => $gff
  });
  # Use variants that all overlap a feature: reference-allele resolution reads
  # the *variant's* slice, which only gets set for variants an annotation source
  # has touched, so a non-overlapping variant would have none.
  my $vcf = $test_cfg->create_input_file([
    ['21', 25585733, 'in_promoter', 'C', 'T'],
    ['21', 25587701, 'in_ctcf',     'T', 'C'],
  ]);
  my $ib = annotate($source, $cfg, $vcf);

  # With only a regulatory source in play nothing has populated
  # transcript_variations, so the output factory would trigger VariationFeature's
  # lazy transcript fetch. A real run always has a transcript source; mark them
  # as computed-and-empty to stand in for one.
  $_->{transcript_variations} ||= {} for @{$ib->buffer};

  # tab
  my $tab = Bio::EnsEMBL::VEP::OutputFactory::Tab->new({config => $cfg});
  my @lines = grep {/RegulatoryFeature/} @{$tab->get_all_lines_by_InputBuffer($ib)};
  ok(scalar @lines, 'tab output - regulatory lines produced');
  like($lines[0], qr/regulatory_region_variant/, 'tab output - consequence');
  like($lines[0], qr/ENSR21_/, 'tab output - feature identifier');
  like($lines[0], qr/RegulatoryFeature/, 'tab output - Feature_type');

  # json
  my $json = Bio::EnsEMBL::VEP::OutputFactory::JSON->new({config => $cfg});
  my $hashes = $json->get_all_output_hashes_by_InputBuffer($ib);
  my ($with_reg) = grep {$_->{regulatory_feature_consequences}} @$hashes;
  ok($with_reg, 'json output - regulatory_feature_consequences present');

  my $rfc = $with_reg->{regulatory_feature_consequences}->[0];
  is($rfc->{regulatory_feature_id}, 'ENSR21_PROM1', 'json output - feature identifier');
  is($rfc->{biotype}, 'promoter', 'json output - biotype carries the feature type');
  ok(
    (grep {$_ eq 'regulatory_region_variant'} @{$rfc->{consequence_terms}}),
    'json output - consequence term'
  );
}


## CONFIG WIRING
################

use_ok('Bio::EnsEMBL::VEP::AnnotationSourceAdaptor');

# --regulatory_gff must turn on the regulatory flag, or Constants.pm drops
# BIOTYPE - the only field distinguishing an enhancer from a promoter
{
  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %{$test_cfg->base_testing_cfg},
    regulatory_gff => $gff,
  });
  is($cfg->param('regulatory'), 1, '--regulatory_gff sets regulatory');
}

{
  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %{$test_cfg->base_testing_cfg},
    regulatory_gff => $gff,
  });
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $cfg});
  my $sources = $asa->get_all_regulatory_gff;
  is(scalar @$sources, 1, 'AnnotationSourceAdaptor creates one source');
  is(
    ref($sources->[0]),
    'Bio::EnsEMBL::VEP::AnnotationSource::File::RegFeat',
    'AnnotationSourceAdaptor creates the right class'
  );
}

{
  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %{$test_cfg->base_testing_cfg},
    regulatory_gff => '/does/not/exist.gff3.gz',
  });
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $cfg});
  throws_ok { $asa->get_all_regulatory_gff } qr/not found/, 'missing file throws';
}

{
  my $cfg = Bio::EnsEMBL::VEP::Config->new({%{$test_cfg->base_testing_cfg}});
  my $asa = Bio::EnsEMBL::VEP::AnnotationSourceAdaptor->new({config => $cfg});
  is(scalar @{$asa->get_all_regulatory_gff}, 0, 'no source without --regulatory_gff');
}

## COEXISTENCE
###############

# Regulatory features must come from exactly one source. merge_features()
# deduplicates within a source only, so a regulatory GFF alongside a regulatory
# cache would double-annotate every overlapping variant.
{
  use_ok('Bio::EnsEMBL::VEP::CacheDir');

  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %{$test_cfg->base_testing_cfg},
    dir            => $test_cfg->{cache_root_dir},
    cache          => 1,
    offline        => 1,
    regulatory_gff => $gff,
  });

  my $cache_dir = Bio::EnsEMBL::VEP::CacheDir->new({
    config => $cfg, dir => $test_cfg->{cache_dir}
  });

  # --regulatory_gff supersedes the cache regulatory source rather than failing.
  # --regulatory_gff itself sets regulatory, so the collision is not one the user
  # asked for; erroring would make the flag unusable with any standard cache,
  # which ships regulatory data.
  my $sources;
  my $warned = warning { $sources = $cache_dir->get_all_AnnotationSources };
  like(
    "$warned", qr/ignoring the regulatory data in this cache/i,
    'regulatory GFF alongside a regulatory cache warns'
  );

  ok(
    !(grep {ref($_) =~ /Cache::RegFeat/} @$sources),
    'cache regulatory source suppressed when --regulatory_gff is set'
  );
  ok(
    (grep {ref($_) =~ /Cache::Transcript/} @$sources),
    'cache transcript source still created alongside a regulatory GFF'
  );
}

# without --regulatory_gff the same cache is fine
{
  my $cfg = Bio::EnsEMBL::VEP::Config->new({
    %{$test_cfg->base_testing_cfg},
    dir        => $test_cfg->{cache_root_dir},
    cache      => 1,
    offline    => 1,
    regulatory => 1,
  });

  my $cache_dir = Bio::EnsEMBL::VEP::CacheDir->new({
    config => $cfg, dir => $test_cfg->{cache_dir}
  });

  my $sources = $cache_dir->get_all_AnnotationSources;
  ok(
    (grep {ref($_) =~ /Cache::RegFeat/} @$sources),
    'regulatory cache source still created without --regulatory_gff'
  );
}


# done
done_testing();
