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
  Bio::EnsEMBL::VEP::AnnotationSource::Cache
  Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(old_maf no_check_alleles failed check_frequency freq_pop freq_freq freq_gt_lt freq_filter)]);

  # checks on frequency filters
  $self->check_frequency_filter if $self->{check_frequency};

  # add this flag to tell VEP runner to use this first
  $self->{can_filter_vfs} = 1;

  return $self;
}

# gets variation cache columns
sub get_cache_columns {
  my $self = shift;

  if(!exists($self->{cols})) {
    my @copy = @Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation::VAR_CACHE_COLS;
    $self->{cols} = \@copy;
    # push @{$self->{cols}}, 'pubmed' if $self->have_pubmed() && $self->param('pubmed');
    # push @{$self->{cols}}, @{$self->{freq_file_pops}} if defined($self->{freq_file_pops});
  }

  return $self->{cols};
}

sub parse_variation {
  my $self = shift;
  my $line = shift;

  my @cols = @{$self->get_cache_columns};
  my $delim = $self->delimiter;
  my @data = split $delim, $line;

  # this switcher is a bit of a hack, should fix cache generation really
  my %v = map {$cols[$_] => $data[$_] eq '.' ? undef : $data[$_]} (0..(@data > @cols ? $#cols : $#data));

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

sub check_frequency_filter {
  my $self = shift;

  throw("ERROR: No population specified") unless $self->{freq_pop};

  my $freq_pop_full = uc($self->{freq_pop});
  $freq_pop_full =~ s/EXAC/ExAC/;
  my ($freq_group, $freq_pop) = split('_', $freq_pop_full);

  my %valid = map {$_ => 1} @{$self->get_valid_populations};

  throw("ERROR: Invalid population ".$self->{freq_pop}." specified") unless $freq_pop_full eq '1KG_ALL' || $valid{$freq_pop_full} || $valid{$freq_pop || ''};
}

sub frequency_check_buffer {
  my $self = shift;
  my $buffer = shift;

  my @passed;
  my $freq_filter = $self->{freq_filter};

  foreach my $vf(@{$buffer->buffer}) {
    $self->get_frequency_data($vf);

    # include only
    if($freq_filter eq 'include') {

      # we can disregard this VF straight away if there's no freq data
      next unless $vf->{_freq_check_pass};

      # all passed?
      if($vf->{_freq_check_all_passed}) {
        push @passed, $vf;
      }
      else {
        # new alt alleles
        my @alts = grep {$vf->{_freq_check_pass}->{$_}} @{$vf->alt_alleles};

        # new allele_string
        $vf->allele_string(join('/', $vf->ref_allele_string, @alts));
      }
    }

    # exclude ones that don't pass
    else {

      # no frequency data so we have to include it
      unless($vf->{_freq_check_pass}) {
        push @passed, $vf;
        next;
      }

      # all failed?
      if($vf->{_freq_check_all_failed}) {
        push @passed, $vf;
      }
      elsif(!$vf->{_freq_check_all_passed}) {
        # new alt alleles
        my @alts = grep {!$vf->{_freq_check_pass}->{$_}} @{$vf->alt_alleles};

        # new allele_string
        $vf->allele_string(join('/', $vf->ref_allele_string, @alts));

        push @passed, $vf;
      }
    }
  }

  my $filtered_count = scalar @{$buffer->buffer} - scalar @passed;

  # replace/update the buffer
  if($filtered_count) {

    # we need to reset so that interval trees etc don't get confused
    $buffer->reset_buffer();

    $buffer->buffer(\@passed);

    # update stats
    $self->stats->increment_filtered_variants($filtered_count);
  }
}

sub _add_check_freq_data_to_vf {
  my $self = shift;
  my $vf = shift;
  my $freq_data = shift;

  my $freq_freq     = $self->{freq_freq};
  my $freq_gt_lt    = $self->{freq_gt_lt};
  my $freq_pop_full = uc($self->{freq_pop});

  my $alts = $vf->alt_alleles();

  my %pass;
  my %used_freqs;
  my $pass_count = 0;

  foreach my $alt(@$alts) {
    my $pass = 0;

    # set freq to 0 if it hasn't been observed
    # this I think is a fair assumption if the position has been fully typed
    my $f = $freq_data->{$alt} || 0;

    $pass = 1 if $f >= $freq_freq and $freq_gt_lt eq 'gt';
    $pass = 1 if $f <= $freq_freq and $freq_gt_lt eq 'lt';
    $used_freqs{$alt} = $f;

    $pass_count += $pass;
    $pass{$alt} = $pass;
  }

  $vf->{_freq_check_freqs}->{$freq_pop_full} = \%used_freqs;
  $vf->{_freq_check_pass} = \%pass;

  # all passed or failed?
  $vf->{_freq_check_all_passed} = 1 if $pass_count == scalar @$alts;
  $vf->{_freq_check_all_failed} = 1 if $pass_count == 0;
}

sub delimiter {
  return qr/ /;
}

sub merge_features {
  return $_[1];
}

1;
