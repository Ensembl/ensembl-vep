=head1 LICENSE

Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::OutputFactory::Tab
#
#

=head1 NAME

Bio::EnsEMBL::VEP::OutputFactory::Tab - VEP format output factory

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::OutputFactory::Tab;

use base qw(Bio::EnsEMBL::VEP::OutputFactory::BaseTab);

use Bio::EnsEMBL::VEP::Utils qw(convert_arrayref);
use Bio::EnsEMBL::VEP::Constants;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    fields
  )]);

  return $self;
}

sub output_hash_to_line {
  my $self = shift;
  my $hash = shift;
  return join("\t", map {defined($hash->{$_}) ? convert_arrayref($hash->{$_}) : '-'} @{$self->fields});
}

sub description_headers {
  my $self = shift;

  my $field_descs = \%Bio::EnsEMBL::VEP::Constants::FIELD_DESCRIPTIONS;

  my %other_descs = map {$_->[0] => $_->[1]} @{$self->get_plugin_headers}, @{$self->get_custom_headers};

  return [
    '## Column descriptions:',
    map {'## '.$_.' : '.($field_descs->{$_} || $other_descs{$_} || '?')}
    @{$self->fields}
  ];
}

sub column_header {
  return '#'.join("\t", @{$_[0]->fields});
}

sub fields {
  my $self = shift;

  if(!defined($self->{fields})) {

    my @fields = (
      @Bio::EnsEMBL::VEP::Constants::DEFAULT_OUTPUT_COLS,
      @{$self->flag_fields},
      map {$_->[0]} (@{$self->get_plugin_headers}, @{$self->get_custom_headers})
    );
    
    $self->{fields} = \@fields;
  }

  return $self->{fields};
}

1;
