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

# EnsEMBL module for Bio::EnsEMBL::VEP::InputBuffer
#
#

=head1 NAME

Bio::EnsEMBL::VEP::InputBuffer - class representing a buffer of VariationFeatures to be processed by VEP

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::InputBuffer;

use parent qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(buffer_size)]);

  my $hashref = $_[0];
  if($hashref) {
    $self->parser($hashref->{parser}) if $hashref->{parser};

    if($hashref->{variation_features}) {
      my $buffer = $self->_full_buffer;
      push @$buffer, @{$hashref->{variation_features}};
    }
  }

  return $self;
}

sub parser {
  my $self = shift;

  if(@_) {
    my $parser = shift;
    assert_ref($parser, 'Bio::EnsEMBL::VEP::Parser');
    $self->{parser} = $parser;
  }

  return $self->{parser};
}

sub next {
  my $self = shift;

  $self->reset_buffer();

  my $full_buffer = $self->_full_buffer();
  my $buffer = $self->buffer();

  my $buffer_size = $self->{buffer_size};

  while(@$full_buffer && @$buffer < $buffer_size) {
    push @$buffer, shift @$full_buffer;
  }

  if(my $parser = $self->parser) {
    while(@$buffer < $buffer_size && (my $vf = $parser->next)) {
      push @$buffer, $vf;
    }
  }

  return $buffer;
}

sub reset_buffer {
  $_[0]->{_buffer} = [];
}

sub buffer {
  return $_[0]->{_buffer} ||= [];
}

sub _full_buffer {
  return $_[0]->{_full_buffer} ||= [];
}

1;