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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation - base class for cache variation sources

=head1 SYNOPSIS

Should not be invoked directly

=head1 DESCRIPTION

Base class for annotation sources that read variant data from a cache.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation;

use Scalar::Util qw(weaken looks_like_number);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Cache
  Bio::EnsEMBL::VEP::AnnotationType::Variation
);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config $config,
                 dir    => string $dir,
               }
  Example    : Not invoked directly
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation
  Exceptions : none
  Caller     : CacheDir
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(old_maf no_check_alleles exclude_null_alleles failed check_frequency freq_pop freq_freq freq_gt_lt freq_filter)]);

  # checks on frequency filters
  $self->check_frequency_filter if $self->{check_frequency};

  # add this flag to tell VEP runner to use this first
  $self->{can_filter_vfs} = 1;

  return $self;
}


=head2 parse_variation

  Arg 1      : string $line_of_data
  Example    : $var = $as->parse_variation($line_of_data);
  Description: Parses a line of data read from the cache into a
               hashref structure representing a known variant.
  Returntype : hashref
  Exceptions : none
  Caller     : read_variations_from_file(), _annotate_cl(), _annotate_pm()
  Status     : Stable

=cut

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

  # # sanity check frequency data
  # foreach my $pop(grep {defined($v{$_})} qw(AFR AMR ASN EAS EUR SAS AA EA)) {
  #   $v{$pop} = undef unless $v{$pop} =~ /^(([ACGTN\-]+?\:)([\-e\.\d]+)\,?)+$/i;
  # }

  return \%v;
}


=head2 get_valid_populations

  Example    : $pops = $as->get_valid_populations();
  Description: Get list of populations for which frequency data is available
               in this cache. Currently a bodge which actually returns all of
               the column names.
  Returntype : hashref
  Exceptions : none
  Caller     : read_variations_from_file(), _annotate_cl(), _annotate_pm()
  Status     : Stable

=cut

sub get_valid_populations {
  return $_[0]->get_cache_columns;
}


=head2 get_frequency_data

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Example    : $as->get_frequency_data($vf);
  Description: Add frequency data read from the cache to the $vf,
               accounting for strand and alleles where necessary.
  Returntype : none
  Exceptions : none
  Caller     : frequency_check_buffer()
  Status     : Stable

=cut

sub get_frequency_data {
  my $self = shift;
  my $vf = shift;

  my $freq_pop_full = uc($self->{freq_pop});
  $freq_pop_full =~ s/EXAC/ExAC/i;
  $freq_pop_full =~ s/gnomad/gnomAD/i;
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
    if(
      $freq_pop_full eq '1KG_ALL' &&
      $ex->{minor_allele} &&
      defined($ex->{minor_allele_freq}) && looks_like_number($ex->{minor_allele_freq})
    ) {
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


=head2 check_frequency_filter

  Example    : $as->check_frequency_filter();
  Description: Checks frequency filter parameters, throws if they
               are invalid.
  Returntype : none
  Exceptions : throws if no or invalid population given
  Caller     : new()
  Status     : Stable

=cut

sub check_frequency_filter {
  my $self = shift;

  throw("ERROR: No population specified") unless $self->{freq_pop};

  my $freq_pop_full = uc($self->{freq_pop});
  $freq_pop_full =~ s/EXAC/ExAC/;
  $freq_pop_full =~ s/ADJ/Adj/;
  $freq_pop_full =~ s/GNOMAD/gnomAD/;
  my ($freq_group, $freq_pop) = split('_', $freq_pop_full);

  my %valid = map {$_ => 1} @{$self->get_valid_populations};

  throw("ERROR: Invalid population ".$self->{freq_pop}." specified")
    unless $freq_pop_full eq '1KG_ALL' ||
    $valid{$freq_pop_full} ||
    ($freq_group eq '1KG' && $valid{$freq_pop || ''});
}


=head2 frequency_check_buffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Example    : $as->frequency_check_buffer($ib);
  Description: Filters the variants in an input buffer by their allele
               frequency according to user parameters and any matching
               known variants. Note that the contents of the InputBuffer
               are modified by this, meaning reset_buffer() is called
               on the InputBuffer object.
  Returntype : none
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

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


=head2 _add_check_freq_data_to_vf

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature $vf
  Arg 2      : hashref $freq_data
  Example    : $as->_add_check_freq_data_to_vf($vf, $freq_data);
  Description: Adds any frequency data used to filter the variant
               to the VariationFeature object using the key
               "_freq_check_freqs"
  Returntype : none
  Exceptions : none
  Caller     : get_frequency_data()
  Status     : Stable

=cut

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

    # set freq to 'NA' if the alternate allele hasn't been observed, or
    # the reference allele is the minor allele and the overlapping variant
    # is multiallelic (more than 1 alternative allele), e.g.:
    # Input: 17:7676154 => G/C
    # Minor allele: G (0.4571)
    # Overlapping variant: rs1042522 (G/C/T)
    my $f = $freq_data->{$alt} || 'NA';

    if ($f eq 'NA') {
      $pass = 0;
    }
    else {
      $pass = 1 if $f >= $freq_freq and $freq_gt_lt eq 'gt';
      $pass = 1 if $f <= $freq_freq and $freq_gt_lt eq 'lt';
    }
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


=head2 delimiter

  Example    : $delim = $as->delimiter();
  Description: Get delimiter used in cache files.
  Returntype : string
  Exceptions : none
  Caller     : parse_variation()
  Status     : Stable

=cut

sub delimiter {
  return qr/ /;
}


=head2 merge_features

  Arg 1      : arrayref $var_hashes
  Example    : $unique_vars = $as->merge_features($var_hashes)
  Description: Compatibility method, returns arg
  Returntype : arrayref
  Exceptions : none
  Caller     : annotate_InputBuffer()
  Status     : Stable

=cut

sub merge_features {
  return $_[1];
}

1;
