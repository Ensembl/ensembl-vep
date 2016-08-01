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

# EnsEMBL module for Bio::EnsEMBL::VEP::Utils
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Utils - VEP utility functions

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Utils;
use Exporter;
use Scalar::Util qw(looks_like_number);

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(
  &format_coords
  &trim_sequences
  &get_time
  &convert_arrayref
  &numberify
  &merge_hashes
);

sub format_coords {
  my ($start, $end) = @_;

  if(defined($start)) {
    if(defined($end)) {
      if($start > $end) {
        return $end.'-'.$start;
      }
      elsif($start == $end) {
        return $start;
      }
      else {
        return $start.'-'.$end;
      }
    }
    else {
      return $start.'-?';
    }
  }
  elsif(defined($end)) {
    return '?-'.$end;
  }
  else  {
    return '-';
  }
}

sub trim_sequences {
  my ($ref, $alt, $start, $end) = @_;

  $start ||= 0;
  $end ||= $start + (length($ref) - 1);

  my $changed = 0;

  # trim from left
  while($ref && $alt && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
    $ref = substr($ref, 1);
    $alt = substr($alt, 1);
    $start++;
    $changed = 1;
  }

  # trim from right
  while($ref && $alt && substr($ref, -1, 1) eq substr($alt, -1, 1)) {
    $ref = substr($ref, 0, length($ref) - 1);
    $alt = substr($alt, 0, length($alt) - 1);
    $end--;
    $changed = 1;
  }

  return [$ref, $alt, $start, $end, $changed];
}

sub convert_arrayref {
  if(ref($_[0]) eq 'ARRAY') {
    return join(($_[1] || ","), @{$_[0]});
  }
  else {
    return $_[0];
  }
}

sub numberify {
  my $ref = shift;
  my $exempt = shift || {};

  if(ref($ref) eq 'HASH') {
    foreach my $k(keys %$ref) {
      if(ref($ref->{$k}) =~ /HASH|ARRAY/) {
        numberify($ref->{$k});
      }
      else {
        $ref->{$k} = $ref->{$k} + 0 if defined($ref->{$k}) && !$exempt->{$k} && looks_like_number($ref->{$k});
      }
    }
  }
  elsif(ref($ref) eq 'ARRAY') {
    foreach my $i(0..((scalar @$ref) - 1)) {
      if(ref($ref->[$i]) =~ /HASH|ARRAY/) {
        numberify($ref->[$i]);
      }
      else {
        $ref->[$i] = $ref->[$i] + 0 if defined($ref->[$i]) && looks_like_number($ref->[$i]);
      }
    }
  }

  return $ref;
}

sub merge_hashes {
  my ($x, $y, $add) = @_;

  foreach my $k (keys %$y) {
    if (!defined($x->{$k})) {
      $x->{$k} = $y->{$k};
    } else {
      if(ref($x->{$k}) eq 'ARRAY') {
        $x->{$k} = merge_arrays($x->{$k}, $y->{$k});
      }
      elsif(ref($x->{$k}) eq 'HASH') {
        $x->{$k} = merge_hashes($x->{$k}, $y->{$k}, $add);
      }
      else {
        $x->{$k} = ($add && looks_like_number($x->{$k}) && looks_like_number($y->{$k}) ? $x->{$k} + $y->{$k} : $y->{$k});
      }
    }
  }
  return $x;
}

sub merge_arrays {
  my ($x, $y) = @_;

  my %tmp = map {$_ => 1} @$x;
  push @$x, grep {!$tmp{$_}} @$y;

  return $x;
}

# gets time
sub get_time {
  my @time = localtime(time());

  # increment the month (Jan = 0)
  $time[4]++;

  # add leading zeroes as required
  for my $i(0..4) {
    $time[$i] = "0".$time[$i] if $time[$i] < 10;
  }

  # put the components together in a string
  my $time =
    ($time[5] + 1900)."-".
    $time[4]."-".
    $time[3]." ".
    $time[2].":".
    $time[1].":".
    $time[0];

  return $time;
}