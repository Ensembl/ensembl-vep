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

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

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
      my $buffer = $self->pre_buffer;
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

  my $pre_buffer = $self->pre_buffer();
  my $buffer = $self->buffer();

  my $buffer_size = $self->{buffer_size};

  while(@$pre_buffer && @$buffer < $buffer_size) {
    push @$buffer, shift @$pre_buffer;
  }

  if(my $parser = $self->parser) {
    while(@$buffer < $buffer_size && (my $vf = $parser->next)) {
      push @$buffer, $vf;
    }
  }

  return $buffer;
}

sub get_cache_regions {
  my $self = shift;
  my $size = shift;

  $size ||= $self->param('cache_region_size');

  if(!exists($self->{temp}->{cache_regions}->{$size})) {
    my @regions = ();
    my %seen = ();

    foreach my $vf(@{$self->buffer}) {
      my $chr = $vf->{chr} || $vf->slice->seq_region_name;
      throw("ERROR: Cannot get chromosome from VariationFeature") unless $chr;

      foreach my $region_start(map {int($vf->{$_} / $size)} qw(start end)) {
        my $key = join(':', ($chr, $region_start));
        next if $seen{$key};

        push @regions, [$chr, $region_start];
        $seen{$key} = 1;
      }
    }

    $self->{temp}->{cache_regions}->{$size} = \@regions;
  }

  return $self->{temp}->{cache_regions}->{$size};
}

sub finish_annotation {
  my $self = shift;
  $_->_finish_annotation for @{$self->buffer};
}

sub reset_buffer {
  $_[0]->{temp} = {};
}

sub buffer {
  return $_[0]->{temp}->{buffer} ||= [];
}

sub pre_buffer {
  return $_[0]->{pre_buffer} ||= [];
}

1;