=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation - database variation annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation;

use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Database
  Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation
);

our @VAR_CACHE_COLS = qw(
  variation_name
  failed
  somatic
  start
  end
  allele_string
  strand
  minor_allele
  minor_allele_freq
  clin_sig
  phenotype_or_disease
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(no_check_alleles failed)]);# check_frequency freq_pop freq_freq freq_gt_lt freq_filter)]);

  # add this flag to tell VEP runner to use this first
  # $self->{can_filter_vfs} = 1;

  # frequency filtering not supported, probably won't be unless someone asks for it
  throw("ERROR: --check_frequency is not supported using database variation annotation source") if $self->param('check_frequency');

  return $self;
}

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;

  my $cache = $self->cache;
  my @return;

  my $cache_region_size = $self->{cache_region_size};

  my $sr_cache = $self->seq_region_cache;

  foreach my $region(@$regions) {

    my ($chr, $region_start) = @{$region};

    my ($s, $e) = (
      ($region_start * $cache_region_size) + 1,
      ($region_start + 1) * $cache_region_size
    );

    # no seq_region_id?
    next unless $sr_cache->{$chr};

    my $phenotype_attrib_id = $self->phenotype_attrib_id || 0;

    my $sth = $self->var_dbc->prepare(qq{
      SELECT
        vf.variation_id, vf.variation_name, IF(fv.variation_id IS NULL, 0, 1),
        vf.somatic, vf.seq_region_start, vf.seq_region_end,
        vf.allele_string, vf.seq_region_strand, vf.minor_allele, vf.minor_allele_freq,
        REPLACE(vf.clinical_significance, " ", "_"),
        IF(FIND_IN_SET(?, evidence_attribs) > 0, 1, 0)
      FROM variation_feature vf
      LEFT JOIN failed_variation fv ON fv.variation_id = vf.variation_id
      WHERE vf.seq_region_id = ?
      AND vf.seq_region_start >= ?
      AND vf.seq_region_start <= ?
    });

    $sth->execute($phenotype_attrib_id, $sr_cache->{$chr}, $s, $e);

    my %v;
    $v{$_} = undef for @VAR_CACHE_COLS;

    my ($var_id, %vars_by_id);
    $sth->bind_col(1, \$var_id);
    $sth->bind_col($_+2, \$v{$VAR_CACHE_COLS[$_]}) for (0..$#VAR_CACHE_COLS);

    my @vars;

    while($sth->fetch) {
      my %v_copy = %v;
      $v_copy{allele_string} =~ s/\s+/\_/g;
      push @vars, \%v_copy if $self->filter_variation(\%v_copy);
    }

    $sth->finish();

    $cache->{$chr}->{$region_start} = \@vars;

    push @return, @vars;
  }

  return \@return;
}

sub seq_region_cache {
  my $self = shift;

  if(!exists($self->{seq_region_cache})) {
    my (%cache, $chr, $id);

    my $sth = $self->var_dbc->prepare(qq{
      SELECT seq_region_id, name FROM seq_region
    });

    $sth->execute();
    $sth->bind_columns(\$id, \$chr);
    $cache{$chr} = $id while $sth->fetch();
    $sth->finish;

    $self->{seq_region_cache} = \%cache;
  }

  return $self->{seq_region_cache};
}

sub have_pubmed {
  my $self = shift;

  if(!defined($self->{have_pubmed})) {
    my $sth = $self->var_dbc->prepare(qq{
      SELECT COUNT(*) FROM variation_citation
    });
    $sth->execute;

    my $count;
    $sth->bind_columns(\$count);
    $sth->fetch();
    $sth->finish();

    $self->{have_pubmed} = $count;
  }

  return $self->{have_pubmed};
}

sub phenotype_attrib_id {
  my $self = shift;

  if(!exists($self->{phenotype_attrib_id})) {
    my $sth = $self->var_dbc->prepare(qq{
      SELECT attrib_id FROM attrib WHERE value = 'Phenotype_or_Disease';
    });

    my $a;
    $sth->execute();
    $sth->bind_columns(\$a);
    $sth->fetch();
    $sth->finish();

    $self->{phenotype_attrib_id} = $a;
  }

  return $self->{phenotype_attrib_id};
}

sub var_dbc {
  my $self = shift;

  if(!exists($self->{var_dbc})) {
    my $va = $self->get_adaptor('variation', 'Variation');
    throw("ERROR: Could not get db object\n") unless $va->db and $va->db->dbc;

    $self->{var_dbc} = $va->db->dbc;
  }

  return $self->{var_dbc};
}

sub merge_features {
  return $_[1];
}

sub info {
  my $self = shift;

  if(!exists($self->{info})) {
    my %info;

    # variation source versions
    if(my $sa = $self->get_adaptor('variation', 'source')) {
      foreach my $source_name(qw(dbSNP COSMIC ClinVar ESP HGMD-PUBLIC)) {
        my $version = $sa->get_source_version($source_name);
        $info{$source_name} = $version if $version;
      }
    }

    $self->{info} = \%info;
  }

  return $self->{info};
}

1;