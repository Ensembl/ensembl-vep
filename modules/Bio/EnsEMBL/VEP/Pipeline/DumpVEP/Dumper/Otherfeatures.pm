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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Otherfeatures;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Core);

sub run {
  my $self = shift;

  my $vep_params = $self->get_vep_params();

  $vep_params->{$_} = 1 for qw(
    protein
    refseq
    all_refseq
    domains
    uniprot
  );

  $vep_params->{'core_type'} = 'otherfeatures';

  my $config = Bio::EnsEMBL::VEP::Config->new($vep_params);

  my $region_size = $self->param('region_size');
  my $hive_dbc = $self->dbc;
  $hive_dbc->disconnect_if_idle() if defined $hive_dbc;

  my $bam = $vep_params->{'bam'} || $self->param('bam');

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Transcript->new({
    config => $config,
    cache_region_size => $region_size,
    bam => $bam
  });

  # bam requires synonyms loaded
  $as->chromosome_synonyms($as->param('synonyms'));

  $config->param($_, 0) for qw(sift polyphen);

  my $cache = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript->new({
    config => $config,
    cache_region_size => $region_size,
    dir => $self->get_cache_dir($vep_params)
  });

  $self->dump_chrs($as, $cache);

  my $extra = {};
  $extra->{bam} = (split('/', $bam))[-1] if $bam;

  $self->dump_info($as, $self->get_cache_dir($vep_params), $extra);
  
  return;
}

1;
