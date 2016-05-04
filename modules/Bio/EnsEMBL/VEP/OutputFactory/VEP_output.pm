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

# EnsEMBL module for Bio::EnsEMBL::VEP::OutputFactory::VEP_output
#
#

=head1 NAME

Bio::EnsEMBL::VEP::OutputFactory::VEP_output - VEP format output factory

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::OutputFactory::VEP_output;

use base qw(Bio::EnsEMBL::VEP::OutputFactory::BaseTab);

use Bio::EnsEMBL::VEP::Utils qw(convert_arrayref);
use Bio::EnsEMBL::VEP::Constants;

my %OUTPUT_COLS_HASH = map {$_ => 1} @Bio::EnsEMBL::VEP::Constants::DEFAULT_OUTPUT_COLS;

sub output_hash_to_line {
  my $self = shift;
  my $hash = shift;

  # "core" fields
  my @line = map {defined($hash->{$_}) ? convert_arrayref($hash->{$_}) : '-'} @Bio::EnsEMBL::VEP::Constants::DEFAULT_OUTPUT_COLS;

  # add additional fields to "Extra" col at the end
  my %extra =
    map {$_ => $hash->{$_}}
    grep {!$OUTPUT_COLS_HASH{$_}}
    keys %$hash;

  my $field_order = $self->field_order;

  push @line, (
    join(';',
      map {$_.'='.convert_arrayref($extra{$_})}
      sort {
        (defined($field_order->{$a}) ? $field_order->{$a} : 100)
        <=>
        (defined($field_order->{$b}) ? $field_order->{$b} : 100)

        ||

        $a cmp $b
      }
      keys %extra
    )
    || '-'
  );

  return join("\t", @line);
}

sub description_headers {
  my $self = shift;

  my $field_descs = \%Bio::EnsEMBL::VEP::Constants::FIELD_DESCRIPTIONS;

  my @headers = '## Column descriptions:';

  push @headers,
    map {'## '.$_.' : '.($field_descs->{$_} || '?')}
    @Bio::EnsEMBL::VEP::Constants::DEFAULT_OUTPUT_COLS;
  
  push @headers, '## Extra column keys:';
  push @headers,
    map {'## '.$_.' : '.($field_descs->{$_} || '?')}
    @{$self->fields};

  return \@headers;
}

sub column_header {
  return '#'.join("\t", (@Bio::EnsEMBL::VEP::Constants::DEFAULT_OUTPUT_COLS, 'Extra'));
}

sub fields {
  my $self = shift;
  return $self->{fields} ||= $self->flag_fields;
}

sub field_order {
  my $self = shift;

  if(!exists($self->{field_order})) {
    my @fields = @{$self->fields};
    
    $self->{field_order}->{$fields[$_]} = $_ for 0..$#fields;
  }

  return $self->{field_order};
}

return 1;
