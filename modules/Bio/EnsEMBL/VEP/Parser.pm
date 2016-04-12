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

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  throw("ERROR: No file given\n") unless $hashref->{file};
  $self->file($hashref->{file});

  $self->line_number(0);

  return $self;
}

sub file {
  my $self = shift;

  if(@_) {
    my $file = shift;

    if(uc($file) eq 'STDIN') {
      $file = uc($file);
    }
    else {
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

1;
