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

# EnsEMBL module for Bio::EnsEMBL::VEP::Utils
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Utils - VEP utility functions

=head1 SYNOPSIS

use Bio::EnsEMBL::VEP::Utils qw(
  format_coords
  convert_arrayref
  merge_hashes
  merge_arrays
);

# 5-10
my $formatted = format_coords(5, 10);

# a,b,c
my $converted = convert_arrayref(['a', 'b', 'c']);

# {
#   a => 1,
#   b => 2,
# }
my $merged = merge_hashes({a => 1}, {b => 2});

# [1, 2, 3]
my $merged2 = merge_arrays([1, 2], [2, 3]);

=head1 DESCRIPTION

The Utils class contains utility methods used by one or more VEP classes
that do not depend on being a class method.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Utils;
use Exporter;
use Scalar::Util qw(looks_like_number);
use FindBin qw($RealBin);
use Bio::EnsEMBL::VEP::Constants;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(
  &format_coords
  &get_time
  &convert_arrayref
  &numberify
  &merge_hashes
  &merge_arrays
  &find_in_ref
  &get_compressed_filehandle
  &get_version_data
  &get_version_string
);

our ($CAN_USE_PERLIO_GZIP, $CAN_USE_GZIP, $CAN_USE_IO_UNCOMPRESS);

BEGIN {

  # check PerlIO::gzip
  if (eval q{ require PerlIO::gzip; 1 }) {
    $CAN_USE_PERLIO_GZIP = 1;
  }

  # check gzip
  if (`which gzip` =~ /\/gzip/) {
    $CAN_USE_GZIP = 1;
  }

  # check IO::Uncompress::Gunzip 
  if(eval q{ use IO::Uncompress::Gunzip qw($GunzipError); 1}) {
    $CAN_USE_IO_UNCOMPRESS = 1;
  }
}


=head2 format_coords

  Arg 1      : (optional) int $start
  Arg 2      : (optional) int $end
  Example    : $formatted = format_coords($start, $end)
  Description: Formats coordinates according to a few rules:
                - $start and $end provided : "$start-$end"
                - $end > $start : "$end-$start"
                - $start == $end : "$start"
                - $start but no $end : "$start-?"
                - $end but no $start : "?-$end"
                - no $start or $end : "?"
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

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


=head2 convert_arrayref

  Arg 1      : arrayref $arrayref or scalar $scalar
  Arg 2      : (optional) string $separator
  Example    : $converted = convert_arrayref($arrayref)
  Description: If given an arrayref, returns string joined by $separator
               (defaults to comma ","). If given a scalar, just returns
               the given scalar.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub convert_arrayref {
  if(ref($_[0]) eq 'ARRAY') {
    return join(($_[1] || ","), @{$_[0]});
  }
  else {
    return $_[0];
  }
}


=head2 numberify

  Arg 1      : arrayref $arrayref or hashref $hashref
  Arg 2      : (optional) hashref $exempt_keys
  Example    : $converted = numberify($hashref)
  Description: For each element of the arrayref or value in the hashref,
               forces any that look like numbers to be stored as an int or float
               rather than a string. Mostly in Perl this makes no difference,
               but it is required when serialising to JSON to properly store
               numbers.

               If the hashref $exempt_keys is provided, any matching keys in the
               input $hashref will not be numberified.
  Returntype : arrayref or hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub numberify {
  my $ref = shift;
  my $exempt = shift || {};

  if(ref($ref) eq 'HASH') {
    foreach my $k(keys %$ref) {
      if(ref($ref->{$k}) =~ /HASH|ARRAY/) {
        numberify($ref->{$k}, $exempt);
      }
      else {
        $ref->{$k} = $ref->{$k} + 0 if defined($ref->{$k}) && !$exempt->{$k} && looks_like_number($ref->{$k});
      }
    }
  }
  elsif(ref($ref) eq 'ARRAY') {
    foreach my $i(0..((scalar @$ref) - 1)) {
      if(ref($ref->[$i]) =~ /HASH|ARRAY/) {
        numberify($ref->[$i], $exempt);
      }
      else {
        $ref->[$i] = $ref->[$i] + 0 if defined($ref->[$i]) && looks_like_number($ref->[$i]);
      }
    }
  }

  return $ref;
}


=head2 merge_hashes

  Arg 1      : hashref $a
  Arg 2      : hashref $b
  Arg 3      : (optional) bool $add_values
  Example    : $merged = merge_hashes($a, $b);
  Description: Iteratively merges hashes $a and $b. By default the value in $a
               is retained if a key is shared between both. If $add_values is
               set to a true value, the values are added together (numerically)
               if both look like a number.
               Uses merge_arrays if any values within the hashrefs are arrayrefs.
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

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


=head2 merge_arrays

  Arg 1      : arrayref $a
  Arg 2      : arrayref $b
  Example    : $merged = merge_arrays($a, $b);
  Description: Merges arrayrefs $a and $b, retaining only unique values.
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub merge_arrays {
  my ($x, $y) = @_;

  my %tmp = map {$_ => 1} @$x;
  push @$x, grep {!$tmp{$_}} @$y;

  return $x;
}


=head2 find_in_ref

  Arg 1      : ref $ref
  Arg 2      : hashref $required_keys
  Arg 3      : (optional) hashref $return
  Arg 4      : (optional) string $this_key
  Example    : $data = find_in_ref($ref, {foo => 1});
  Description: Find values for the keys specified in $required_keys in the
               arbitrarily nested data structure $ref. Values found for
               the required keys are added to the $return hashref as:
               {
                 key1 => [value1],
                 key2 => [value2, value3]
               }

               Optional args $return and $this_key are used internally as
               the method runs recursively.
  Returntype : hashref
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub find_in_ref {
  my ($ref, $want_keys, $return, $this_key) = @_;

  $return ||= {};

  if(ref($ref) eq 'HASH') {
    foreach my $key(keys %$ref) {
      if(ref($ref->{$key})) {
        find_in_ref($ref->{$key}, $want_keys, $return, $key);
      }
      else {
        merge_arrays($return->{$key} ||= [], [$ref->{$key}]) if $want_keys->{$key};
      }
    }
  }
  elsif(ref($ref) eq 'ARRAY') {
    foreach my $el(@$ref) {
      if(ref($el)) {
        find_in_ref($el, $want_keys, $return);
      }
      else {
        merge_arrays($return->{$this_key} ||= [], [$el]) if $this_key && $want_keys->{$this_key};
      }
    } 
  }
  elsif($this_key && $want_keys->{$this_key}) {
    merge_arrays($return->{$this_key} ||= [], [$ref]);
  }

  return $return;
}


=head2 get_compressed_filehandle

  Arg 1      : string $filename
  Arg 2      : (optional) bool $multistream
  Example    : $fh = get_compressed_filehandle($file);
  Description: Gets a filehandle for reading from a compressed input file.
               Uses one of the following in order of preference:
               PerlIO::gzip, gzip, IO::Uncompress::Gunzip

               If $multistream is set to true, then PerlIO::gzip is skipped.
               $multistream must be set for block-gzipped files e.g. those
               compressed with bgzip.
  Returntype : glob
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_compressed_filehandle {
  my ($file, $multi) = @_;

  die("ERROR: No file given\n") unless $file;
  die("ERROR: File $file does not exist\n") unless -e $file;
  die("ERROR: File $file does not look like a binary file\n") unless -B $file;

  if($CAN_USE_PERLIO_GZIP && !$multi) {
    open my $fh, "<:gzip", $file or die("ERROR: $!");
    return $fh;
  }
  elsif($CAN_USE_GZIP) {
    open my $fh, "gzip -dc $file |" or die("ERROR: $!");
    return $fh;
  }
  elsif($CAN_USE_IO_UNCOMPRESS) {
    my $fh;
    eval q{ $fh = IO::Uncompress::Gunzip->new($file, MultiStream => $multi) or die("ERROR: $GunzipError"); };
    die($@) if $@;
    return $fh;
  }
  else {
    die("Cannot read from compressed or binary file");
  }
}


=head2 get_time

  Example    : $time = get_time();
  Description: Gets a timestamp of the current time.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

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


=head2 get_version_data

  Arg 1      : (optional) string $dir_with_version_data
  Example    : $version_data = get_version_data($dir)
  Description: Gets a hashref of version data for VEP and any
               Ensembl API modules installed by VEP's INSTALL.pl
               script. $dir_with_version_data defaults to
               "ensembl-vep/.version/"
  Returntype : hashref
  Exceptions : none
  Caller     : get_version_string(), INSTALL.pl
  Status     : Stable

=cut

sub get_version_data {
  my $version_dir = shift || $RealBin.'/.version';
  my %version_data = ();

  $version_data{'ensembl-vep'} = {
    'release' => $Bio::EnsEMBL::VEP::Constants::VEP_VERSION,
    'sub'     => $Bio::EnsEMBL::VEP::Constants::VEP_SUB_VERSION,
  };

  opendir DIR, $version_dir or return \%version_data;
  foreach my $module(grep {!/^\./} readdir(DIR)) {
    my %hash = ();

    open IN, $version_dir.'/'.$module;
    while(<IN>) {
      chomp;
      my @data = split;
      $hash{$data[0]} = $data[1];
    }
    close IN;

    $version_data{$module} = \%hash;
  }
  closedir DIR;

  return \%version_data;
}


=head2 get_version_string

  Arg 1      : (optional) string $dir_with_version_data
  Example    : $version_string = get_version_string($dir)
  Description: Gets a string of version data for VEP and any
               Ensembl API modules installed by VEP's INSTALL.pl
               script, suitable for printing as part of a help
               message. $dir_with_version_data defaults to
               "ensembl-vep/.version/"
  Returntype : string
  Exceptions : none
  Caller     : vep, haplo
  Status     : Stable

=cut

sub get_version_string {
  my $version_data = get_version_data(@_);
  return join(
    "\n  ",
    map {
      sprintf("%-20s : %s", $_, $version_data->{$_}->{release}).
      (defined($version_data->{$_}->{sub}) ? '.'.substr($version_data->{$_}->{sub}, 0, 7) : '')
    }
    sort keys %$version_data
  );
}

1;

