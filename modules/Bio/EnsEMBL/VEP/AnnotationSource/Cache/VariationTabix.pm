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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix - class for cache variation source indexed with tabix

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::Cache::VariationTabix;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(
  Bio::EnsEMBL::VEP::AnnotationSource::Cache::BaseCacheVariation
  Bio::EnsEMBL::VEP::AnnotationSource::BaseVariation
);

our $CAN_USE_TABIX_PM;

BEGIN {
  if (eval { require Bio::DB::HTS::Tabix; 1 }) {
    $CAN_USE_TABIX_PM = 1;
  }
  else {
    $CAN_USE_TABIX_PM = 0;
  }
}

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

# uses command line tabix util
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

      my $region_string = join " ", map {
        $source_chr.':'.($_->{start} > $_->{end} ?
        $_->{end}.'-'.$_->{start} :
        $_->{start}.'-'.$_->{end})
      } @tmp_list;

      open VARS, "tabix -f $file $region_string 2>&1 |"
        or die "\nERROR: Could not open tabix pipe for $file\n";

      # convert list to hash so we can look up quickly by position
      my %hash;
      push @{$hash{$_->{start}}}, $_ for @tmp_list;

      VAR: while(<VARS>) {
        chomp;
        my $existing = $self->parse_variation($_);

        foreach my $input(@{$hash{$existing->{start}} || []}) {
          if(
            $self->filter_variation($existing) &&
            !$self->is_var_novel($existing, $input)
          ) {
            push @{$input->{existing}}, $existing unless
              grep {$_->{variation_name} eq $existing->{variation_name}}
              @{$input->{existing} || []};
          }
        }
      }

      close VARS;

      $_->{existing} ||= [] for @tmp_list;
    }
  }
}

# uses tabix perl module
sub _annotate_pm {
  my $self = shift;
  my $by_chr = shift;

  my $p = 0;

  foreach my $chr(keys %$by_chr) {

    my $source_chr = $self->get_source_chr_name($chr);
    my $file = $self->get_dump_file_name($source_chr);
    next unless -e $file;
    my $tabix_obj = $self->{tabix_obj}->{$chr} ||= Bio::DB::HTS::Tabix->new(filename => $file);
    next unless $tabix_obj;

    foreach my $vf(@{$by_chr->{$chr}}) {
      my $iter = $tabix_obj->query(sprintf("%s:%i-%i", $source_chr, $vf->{start} - 1, $vf->{end} + 1));
      next unless $iter;

      while(my $line = $iter->next) {
        chomp $line;
        my $existing = $self->parse_variation($line);

        if(
          $self->filter_variation($existing) &&
          !$self->is_var_novel($existing, $vf)
        ) {
          push @{$vf->{existing}}, $existing unless
            grep {$_->{variation_name} eq $existing->{variation_name}}
            @{$vf->{existing} || []};
        }
      }
    }
  }
}

sub delimiter {
  return "\t";
}

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
