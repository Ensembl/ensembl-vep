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

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(
  &format_coords
  &get_time
  &convert_arrayref
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

sub convert_arrayref {
  if(ref($_[0]) eq 'ARRAY') {
    return join(($_[1] || ","), @{$_[0]});
  }
  else {
    return $_[0];
  }
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