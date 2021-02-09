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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper;

use strict;
use warnings;

use Storable qw(nstore_fd);
use File::Path qw(mkpath);

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::BaseVEP);

sub param_defaults {
  return {
    'sift'     => 0,
    'polyphen' => 0,
    'convert'  => 0,
    'is_multispecies' => 0,
    'dir_suffix' => '',
    'eg_version' => undef,
  };
}

sub get_vep_params {
  my $self = shift;

  my $params = {};

  # basic params
  $params->{dir}     = $self->data_dir;
  $params->{debug}   = $self->param('debug');
  $params->{species} = $self->required_param('species');
  $params->{host}    = $self->required_param('host');
  $params->{port}    = $self->required_param('port');
  $params->{user}    = $self->required_param('user');
  $params->{pass}    = $self->required_param('pass') || '';

  $params->{db_version}      = $self->required_param('db_version');
  $params->{cache_version}   = $self->param('eg_version') || $self->param('ensembl_release');
  $params->{assembly}        = $self->required_param('assembly');
  $params->{is_multispecies} = $self->param('is_multispecies');
  $params->{no_slice_cache}  = 1;

  # add synonyms file, needed by BAM
  $params->{synonyms} = sprintf(
    '%s/synonyms/%s_%s_chr_synonyms.txt',
    $self->param('pipeline_dump_dir'),
    $self->param('species'),
    $self->param('assembly')
  );

  # sift, polyphen
  $params->{$_} = 'b' for grep {$self->param($_)} qw(sift polyphen);

  # species-specific
  my $species_flags = $self->param('species_flags');
  
  if(my $flags = $species_flags->{$params->{species}}) {
    
    # assembly-specific
    if(my $as = $flags->{assembly_specific}) {
      delete $flags->{assembly_specific};

      my $assembly = $params->{assembly};
      
      if(my $as_flags = $as->{$assembly}) {
        $params->{$_} = $as_flags->{$_} for keys %$as_flags;
      }
    }

    $params->{$_} = $flags->{$_} for keys %$flags;
  }

  return $params;
}

sub get_cache_dir {
  my ($self, $vep_params) = @_;

  return sprintf(
    '%s/%s/%s_%s',
    $vep_params->{dir},
    $vep_params->{species}.($vep_params->{refseq} ? '_refseq' : ''),
    $vep_params->{cache_version},
    $vep_params->{assembly}
  );
}

sub dump_chrs {
  my ($self, $db_as, $cache_as) = @_;

  my $sa = $db_as->get_adaptor('core', 'slice');

  my @regions = @{$self->param('regions')};

  my $region_size = $self->param('region_size');

  while(my $region = shift @regions) {
    my ($chr, $sr, $slice_start, $slice_end) = (
      $region->{chr},
      $region->{seq_region_id},
      $region->{start},
      $region->{end}
    );

    my $slice = $sa->fetch_by_seq_region_id($sr);

    my $s = int($slice_start / $region_size);
    my $l = int($slice_end / $region_size);
    my $first = 1;

    while($s <= $l) {
      # if($s != 89) { $s++; next; }
      my $obj = $self->get_dumpable_object($db_as, $sr, $chr, $s);
      my $file = $cache_as->get_dump_file_name($chr, ($s  * $region_size) + 1, ($s + 1) * $region_size);

      if($first) {
        my $filedir = $file;
        $filedir =~ s/\/[^\/]+$//;
        mkpath($filedir) unless -d $filedir;
        $first = 0;
      }

      $self->dump_obj($obj, $file, $chr);

      $self->post_dump($obj, $db_as, $chr);

      $db_as->clean_cache();
      $s++;
    }
  }
}

sub dump_obj {
  my $self = shift;
  my $obj = shift;
  my $file = shift;

  open my $fh, "| gzip -9 -c > ".$file or die "ERROR: Could not write to dump file $file";
  nstore_fd($obj, $fh);
  close $fh;
}

sub post_dump {
  return 1;
}

1;
