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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation - database variation annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
  config => $config,
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Database-based annotation source for known variant data.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation;

use Scalar::Util qw(weaken);
use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Database
  Bio::EnsEMBL::VEP::AnnotationType::Variation
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


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config $config,
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation
  Exceptions : throws if --check_frequency set (only compatible with cache annotation sources)
  Caller     : AnnotationSourceAdaptor
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(no_check_alleles exclude_null_alleles failed)]);# check_frequency freq_pop freq_freq freq_gt_lt freq_filter)]);

  # add this flag to tell VEP runner to use this first
  # $self->{can_filter_vfs} = 1;

  # frequency filtering not supported, probably won't be unless someone asks for it
  throw("ERROR: --check_frequency is not supported using database variation annotation source") if $self->param('check_frequency');

  return $self;
}


=head2 get_features_by_regions_uncached

  Arg 1      : arrayref $regions
  Example    : $trs = $as->get_features_by_regions_uncached($regions)
  Description: Gets all known variants overlapping the given set of regions. See
               Bio::EnsEMBL::VEP::AnnotationSource::get_all_regions_by_InputBuffer()
               for information about regions.
  Returntype : arrayref of variant hashrefs
  Exceptions : none
  Caller     : get_all_features_by_InputBuffer()
  Status     : Stable

=cut

sub get_features_by_regions_uncached {
  my $self = shift;
  my $regions = shift;
  my $chr_is_seq_region = shift;

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

    my $adaptor = $self->get_adaptor('variation', 'phenotypefeature');
    my $source_id = $self->clinvar_source_id_cache;
    my $attribs = $adaptor->get_clinsig_alleles_by_location($chr_is_seq_region ? $chr : $sr_cache->{$chr}, $s, $e, $source_id) if defined($adaptor) && defined($source_id);

    # no seq_region_id?
    next unless $sr_cache->{$chr} || $chr_is_seq_region;

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
    }, {mysql_use_result => 1});

    $sth->execute($phenotype_attrib_id, $chr_is_seq_region ? $chr : $sr_cache->{$chr}, $s, $e);

    my %v;
    $v{$_} = undef for @VAR_CACHE_COLS;

    my ($var_id, %vars_by_id);
    $sth->bind_col(1, \$var_id);
    $sth->bind_col($_+2, \$v{$VAR_CACHE_COLS[$_]}) for (0..$#VAR_CACHE_COLS);

    my @vars;
    while($sth->fetch) {
      my %v_copy = %v;
      $v_copy{allele_string} =~ s/\s+/\_/g;
      my $v_clinsigs = $attribs->{($chr_is_seq_region ? $chr : $sr_cache->{$chr}) . ':' . $v_copy{start} . '-' . $v_copy{end}};
      my @pfas_by_allele;
      my %clin_sigs;
      foreach my $pfa(@{$v_clinsigs})
      {
	if(defined($pfa->{clinvar_clin_sig}) && $v_copy{variation_name} eq $pfa->{id})
	{
          $pfa->{clinvar_clin_sig}=~s/ /_/g;
          $clin_sigs{$pfa->{risk_allele} . ':' .$pfa->{clinvar_clin_sig}} = 1;
        }
      }
      my @array = keys(%clin_sigs);
      $v_copy{clin_sig_allele} = join ';', @array if scalar(@array);
      $v_copy{variation_id} = $var_id;
      ## fix for e!94 alleles
      $v_copy{allele_string} =~ s/\/$//g;
      $v_copy{allele_string} =~ s/\/\//\//;

      push @vars, \%v_copy if $self->filter_variation(\%v_copy);
    }

    $sth->finish();

    $cache->{$chr}->{$region_start} = \@vars;

    push @return, @vars;
  }

  return \@return;
}


=head2 seq_region_cache

  Example    : $cache = $as->seq_region_cache()
  Description: Gets a hashref mapping chromosome name to the database's
               internal seq_region_id.
  Returntype : hashref
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

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

=head2 clinvar_source_id_cache

  Example    : $id = $as->clinvar_source_id_cache()
  Description: Gets the id for the database's
               internal source_id.
  Returntype : scalar
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub clinvar_source_id_cache {
  my $self = shift;

  if(!exists($self->{clinvar_source_id_cache})) {
    my ($id);

    my $sth = $self->var_dbc->prepare(qq{
      select source_id from source where name = 'ClinVar'
    });

    $sth->execute();
    $sth->bind_columns(\$id);
    $self->{clinvar_source_id_cache} = $id  while $sth->fetch();
    $sth->finish;

  }

  return $self->{clinvar_source_id_cache};
}






=head2 have_pubmed

  Example    : $have_pubmed = $as->have_pubmed()
  Description: Finds if this database contains any database citation information.
  Returntype : bool
  Exceptions : none
  Caller     : DumpVEP pipeline
  Status     : Stable

=cut

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


=head2 phenotype_attrib_id

  Example    : $id = $as->phenotype_attrib_id()
  Description: Gets the internal database attribute ID that indicates if
               a variant has an associated phenotype or disease.
  Returntype : int
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

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


=head2 phenotype_attrib_id

  Example    : $dbc = $as->var_dbc()
  Description: Gets the DBI database connection object for the variation database.
  Returntype : DBI
  Exceptions : none
  Caller     : get_features_by_regions_uncached()
  Status     : Stable

=cut

sub var_dbc {
  my $self = shift;

  if(!exists($self->{var_dbc})) {
    my $va = $self->get_adaptor('variation', 'Variation');
    throw("ERROR: Could not get db object\n") unless $va->db and $va->db->dbc;

    $self->{var_dbc} = $va->db->dbc;
  }

  return $self->{var_dbc};
}


=head2 merge_features

  Arg 1      : arrayref $variant_hashrefs
  Example    : $variant_hashrefs = $as->merge_features($variant_hashrefs)
  Description: Compatibility method, returns arg
  Returntype : arrayref
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub merge_features {
  return $_[1];
}


=head2 info

  Example    : $info = $as->info()
  Description: Gets the info hashref for this annotation source. Contains
               version information for:
                - dbSNP
                - COSMIC
                - ClinVar
                - ESP
                - HGMD-PUBLIC
  Returntype : hashref
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::BaseRunner
  Status     : Stable

=cut

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
