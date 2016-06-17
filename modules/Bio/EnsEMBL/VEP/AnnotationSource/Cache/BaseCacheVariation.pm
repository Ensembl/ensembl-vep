=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation - base class for cache variation sources

DO NOT USE DIRECTLY

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation
  Bio::EnsEMBL::VEP::AnnotationSource::Cache
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(old_maf)]);

  # checks on frequency filters
  $self->check_frequency_filter if $self->{check_frequency};

  return $self;
}

sub parse_variation {
  my $self = shift;
  my $line = shift;

  my @cols = @{$self->get_cache_columns};
  my $delim = $self->delimiter;
  my @data = split $delim, $line;

  # assumption fix for old cache files
  push @cols, ('AFR', 'AMR', 'ASN', 'EUR') if scalar @data > scalar @cols;

  my %v = map {$cols[$_] => $data[$_] eq '.' ? undef : $data[$_]} (0..$#data);

  $v{$_} ||= 0 for qw(failed somatic phenotype_or_disease);
  $v{end}     ||= $v{start};
  $v{strand}  ||= 1;

  # hack for odd frequency data
  if($self->{old_maf}) {
    foreach my $pop(grep {defined($v{$_})} qw(AFR AMR ASN EUR)) {
     $v{$pop} =~ s/^.+?\://;
     $v{$pop} =~ s/\,.+//g;
     $v{$pop} = 1 - $v{$pop} if $v{$pop} =~ /\d+/ && $v{$pop} > 0.5;
    }
  }

  # sanity check frequency data
  foreach my $pop(grep {defined($v{$_})} qw(AFR AMR ASN EAS EUR SAS AA EA)) {
    $v{$pop} = undef unless $v{$pop} =~ /^([ACGTN-]+\:)?(0|0\.\d+|1)$/;
  }

  return \%v;
}

sub get_valid_populations {
  return $_[0]->get_cache_columns;
}

sub get_frequency_data {
  my $self = shift;
  my $vf = shift;

  my $freq_pop_full = uc($self->{freq_pop});
  $freq_pop_full =~ s/EXAC/ExAC/;
  my ($freq_group, $freq_pop) = split('_', $freq_pop_full);

  my $vf_strand = $vf->{strand} || 1;

  my %freq_data;

  foreach my $ex(@{$vf->{existing}}) {
    my $ex_strand = $ex->{strand} || 1;

    my %ex_alleles;
    foreach my $a(split('/', $ex->{allele_string})) {
      reverse_comp(\$a) unless $vf_strand == $ex_strand;
      $ex_alleles{$a} = 1;
    }

    my $total_freq = 0;

    # special case 1KG_ALL
    if($freq_pop_full eq '1KG_ALL' && defined($ex->{minor_allele_freq}) && defined($ex->{minor_allele})) {
      my $a = $ex->{minor_allele};
      my $f = $ex->{minor_allele_freq};
      
      reverse_comp(\$a) unless $vf_strand == $ex_strand;
      $freq_data{$a} = $f;

      # track total and delete observed alleles
      $total_freq += $f;
      delete($ex_alleles{$a});
    }
    # otherwise check match on first full pop name (e.g. ExAC_AMR, then last part e.g. AMR)
    elsif(my $tmp = ($ex->{$freq_pop_full} || $ex->{$freq_pop || ''})) {
      foreach my $a_f(split(',', $tmp)) {
        my ($a, $f) = split(':', $a_f);

        reverse_comp(\$a) unless $vf_strand == $ex_strand;
        $freq_data{$a} = $f;
        
        # track total and delete observed alleles
        $total_freq += $f;
        delete($ex_alleles{$a});
      }
    }

    # "fill in" by subtracting
    # but we can only do this if we have one remaining allele
    # (typically the REF in a biallelic SNP)
    if(scalar keys %ex_alleles == 1) {
      my $a = (keys %ex_alleles)[0];
      $freq_data{$a} = 1 - $total_freq;
    }
  }

  $self->_add_check_freq_data_to_vf($vf, \%freq_data) if %freq_data;
}

sub delimiter {
  return qr/ /;
}

sub merge_features {
  return $_[1];
}

1;
