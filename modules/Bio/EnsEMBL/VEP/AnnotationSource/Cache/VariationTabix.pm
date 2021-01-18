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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix - class for cache variation source indexed with tabix

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix->new({
  config => $config,
  dir    => $dir
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Cache-based (tabix-converted) annotation source for known variant data.

Variants are stored in one file per chromosome, named [chr]/all_vars.gz

Similar format as that used by Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation,
with some changes:

 - tab "\t" used as delimiter (required by tabix)
 - empty fields replaced by "."
 - chromosome column added to beginning (required by tabix)

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix;

use Scalar::Util qw(weaken);
use FindBin qw($RealBin);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation
);

our ($CAN_USE_TABIX_PM, $CAN_USE_TABIX_CL, $TABIX_BIN);

BEGIN {
  if (eval q{ require Bio::DB::HTS::Tabix; 1 }) {
    $CAN_USE_TABIX_PM = 1;
  }

  $TABIX_BIN = `which tabix`;
  chomp($TABIX_BIN);
  $TABIX_BIN ||= "$RealBin/htslib/tabix";

  if (-e $TABIX_BIN) {
    $CAN_USE_TABIX_CL = 1;
  }
}


=head2 annotate_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer
  Example    : $as->annotate_InputBuffer($ib);
  Description: Gets overlapping known variants for the variants in
               the input buffer, checks if they are "novel" (see is_var_novel())
               and adds them to the key "existing" on the VariationFeature
               object. May also filter the input buffer by comparing frequencies
               of known variants to user-specified thresholds.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub annotate_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  # we only care about non-SVs here
  my %by_chr;
  push @{$by_chr{$_->{chr}}}, $_ for grep {ref($_) eq 'Bio::EnsEMBL::Variation::VariationFeature'} @{$buffer->buffer};

  if($CAN_USE_TABIX_PM) {
    $self->_annotate_pm(\%by_chr);
  }
  else {
    $self->_annotate_cl(\%by_chr);
  }

  $self->frequency_check_buffer($buffer) if $self->{check_frequency};
}


=head2 _annotate_cl

  Arg 1      : hashref $vfs_by_chr
  Example    : $as->_annotate_cl($vfs_by_chr);
  Description: Uses command-line tabix tool to read variants from
               cache files and adds known variants to input variants
               under the "existing" key.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub _annotate_cl {
  my $self = shift;
  my $by_chr = shift;

  my $max = 200;
  my $p = 0;

  foreach my $chr(keys %$by_chr) {
    my $list = $by_chr->{$chr};
    my $source_chr = $self->get_source_chr_name($chr);
    my $file = $self->get_dump_file_name($source_chr);
    next unless -e $file;

    while(scalar @$list) {
      my @tmp_list = sort {$a->{start} <=> $b->{start}} splice @$list, 0, $max;
      $p += scalar @tmp_list;

      my $max_length = 1;
      my @regions;

      foreach my $var(@tmp_list) {
        my $l = ($var->{end} - $var->{start}) + 1;
        $max_length = $l if $l > $max_length;
        push @regions, $source_chr.':'.($var->{start} - 1).'-'.($var->{end} + 1);
      }
      my $region_string = join " ", @regions;

      open VARS, "$TABIX_BIN -f $file $region_string 2>&1 |"
        or throw "\nERROR: Could not open tabix pipe for $file\n";

      # convert list to hash so we can look up quickly by position
      my %hash;
      push @{$hash{$_->{start}}}, $_ for @tmp_list;

      VAR: while(<VARS>) {
        chomp;
        my $existing = $self->parse_variation($_);
        next unless $self->filter_variation($existing);

        for my $start(grep {defined($hash{$_})} (($existing->{start} - $max_length)..($existing->{end} + $max_length))) {
          foreach my $vf(@{$hash{$start}}) {
            next if grep {$_->{variation_name} eq $existing->{variation_name}} @{$vf->{existing} || []};

            my $matched = $self->compare_existing($vf, $existing);
            push @{$vf->{existing}}, $matched if $matched;
          }
        }
      }

      close VARS;

      $_->{existing} ||= [] for @tmp_list;
    }
  }
}

=head2 _annotate_pm

  Arg 1      : hashref $vfs_by_chr
  Example    : $as->_annotate_pm($vfs_by_chr);
  Description: Uses Bio::DB::HTS::Tabix to read variants from
               cache files and adds known variants to input variants
               under the "existing" key. Preferred to _annotate_cl()
               for speed.
  Returntype : none
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub _annotate_pm {
  my $self = shift;
  my $by_chr = shift;

  my $p = 0;

  foreach my $chr(keys %$by_chr) {

    my $source_chr = $self->get_source_chr_name($chr);
    my $tabix_obj = $self->_get_tabix_obj($source_chr);
    next unless $tabix_obj;

    foreach my $vf(@{$by_chr->{$chr}}) {
      my $iter = $tabix_obj->query(sprintf("%s:%i-%i", $source_chr, $vf->{start} - 1, $vf->{end} + 1));
      next unless $iter;

      while(my $line = $iter->next) {
        chomp $line;
        my $existing = $self->parse_variation($line);
        next unless $self->filter_variation($existing);

        next if grep {$_->{variation_name} eq $existing->{variation_name}} @{$vf->{existing} || []};

        my $matched = $self->compare_existing($vf, $existing);
        push @{$vf->{existing}}, $matched if $matched;
      }
    }
  }
}


=head2 _get_tabix_obj

  Arg 1      : string $chr
  Example    : $as->_get_tabix_obj($chr);
  Description: Get Bio::DB::HTS::Tabix object for this chromosome.
               Uses a cache that limits the number of open filehandles.
  Returntype : Bio::DB::HTS::Tabix
  Exceptions : none
  Caller     : _annotate_pm()
  Status     : Stable

=cut

sub _get_tabix_obj {
  my ($self, $chr) = @_;

  # use a cache and limit the number of open files
  my $cache = $self->{_tabix_obj_cache} ||= [];
  my $tabix_obj;

  unless(($tabix_obj) = map {$_->{obj}} grep {$_->{chr} eq $chr} @$cache) {
    my $file = $self->get_dump_file_name($chr);
    
    if($file && -e $file) {
      $tabix_obj = Bio::DB::HTS::Tabix->new(filename => $file);
    }

    push @$cache, { obj => $tabix_obj, chr => $chr };

    # restrict number of open objects
    while(scalar @$cache > 5) {
      my $tmp_hash = shift @$cache;
      $tmp_hash->{obj}->close() if $tmp_hash->{obj};
    }
  }

  return $tabix_obj;
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
  return "\t";
}


=head2 get_dump_file_name

  Arg 1      : string $chr
  Example    : $file = $as->get_dump_file_name(1);
  Description: Gets file name from the cache given a chromosome name.
  Returntype : string
  Exceptions : none
  Caller     : _annotate_cl(), _annotate_pm()
  Status     : Stable

=cut

sub get_dump_file_name {
  my $self = shift;
  my $chr  = shift;

  throw("No chromosome given") unless $chr;

  return sprintf(
    "%s/%s/all_vars\.gz",
    $self->dir,
    $chr,
  );
}

1;
