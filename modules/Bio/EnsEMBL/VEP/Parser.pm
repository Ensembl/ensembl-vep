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

# EnsEMBL module for Bio::EnsEMBL::VEP::Parser
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Parser - base class for input parsers

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Parser;

use parent qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::VEP::Parser::VCF;
use Bio::EnsEMBL::VEP::Parser::VEP_input;
use Bio::EnsEMBL::VEP::Parser::ID;
use Bio::EnsEMBL::VEP::Parser::HGVS;

use Scalar::Util qw(openhandle);
use FileHandle;

my %FORMAT_MAP = (
  'vcf'     => 'VCF',
  'ensembl' => 'VEP_input',
  'id'      => 'ID',
  'hgvs'    => 'HGVS',
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  throw("ERROR: No file given\n") unless $hashref->{file};
  $self->file($hashref->{file});

  $self->line_number(0);

  if(my $format = $hashref->{format}) {

    delete $hashref->{format};

    # detect format
    if(lc($format) eq 'guess' || lc($format) eq 'detect' || lc($format) eq 'auto') {
      $format = $self->detect_format();
    }

    throw("ERROR: Can't detect format\n") unless $format;

    $format = lc($format);
    throw("ERROR: Unknown or unsupported format $format\n") unless $FORMAT_MAP{$format};

    my $class = 'Bio::EnsEMBL::VEP::Parser::'.$FORMAT_MAP{$format};
    return $class->new({%$hashref, config => $self->config});
  }

  return $self;
}

sub file {
  my $self = shift;

  if(@_) {
    my $file = shift;

    if(uc($file) eq 'STDIN') {
      $file = *STDIN;
    }
    elsif(!openhandle($file)) {
      throw("ERROR: File \"$file\" does not exist\n") unless -e $file;
    }

    $self->{file} = $file;
  }
  
  return $self->{file};
}

sub line_number {
  my $self = shift;
  $self->{line_number} = shift if @_;
  return $self->{line_number};
}

sub detect_format {
  my $self = shift;

  my $file = $self->file;

  my $fh;
  my $file_was_fh = 0;

  if(openhandle($file)) {
    $fh = $file;
    $file_was_fh = 1;
  }
  else {
    $fh = FileHandle->new();
    $fh->open($file) or throw("ERROR: Could not open $file to detect format\n");
  }

  my $format;

  while(<$fh>) {
    next if /^\#/;
    chomp;

    my @data = split /\s+/, $_;

    # HGVS: ENST00000285667.3:c.1047_1048insC
    if (
      scalar @data == 1 &&
      $data[0] =~ /^([^\:]+)\:.*?([cgmrp]?)\.?([\*\-0-9]+.*)$/i
    ) {
      $format = 'hgvs';
    }

    # variant identifier: rs123456
    elsif (
      scalar @data == 1
    ) {
      $format = 'id';
    }

    # VCF: 20  14370  rs6054257  G  A  29  0  NS=58;DP=258;AF=0.786;DB;H2  GT:GQ:DP:HQ
    elsif (
      $data[0] =~ /(chr)?\w+/ &&
      $data[1] =~ /^\d+$/ &&
      $data[3] && $data[3] =~ /^[ACGTN\-\.]+$/i &&
      $data[4]
    ) {

      # do some more thorough checking on the ALTs
      my $ok = 1;

      foreach my $alt(split(',', $data[4])) {
        $ok = 0 unless $alt =~ /^[\.ACGTN\-\*]+$|^(\<[\w\:\*]+\>)$/i;
      }

      $format = 'vcf' if $ok;
    }

    # pileup: chr1  60  T  A
    elsif (
      $data[0] =~ /(chr)?\w+/ &&
      $data[1] =~ /^\d+$/ &&
      $data[2] && $data[2] =~ /^[\*ACGTN-]+$/i &&
      $data[3] && $data[3] =~ /^[\*ACGTNRYSWKM\+\/-]+$/i
    ) {
      $format = 'pileup';
    }

    # ensembl: 20  14370  14370  A/G  +
    elsif (
      $data[0] =~ /\w+/ &&
      $data[1] =~ /^\d+$/ &&
      $data[2] && $data[2] =~ /^\d+$/ &&
      $data[3] && $data[3] =~ /(ins|dup|del)|([ACGTN-]+\/[ACGTN-]+)/i
    ) {
      $format = 'ensembl';
    }

    # reset file handle if it was a handle
    seek $fh, 0, 0 if $file_was_fh;
    last;
  }

  return $format;
}

1;
