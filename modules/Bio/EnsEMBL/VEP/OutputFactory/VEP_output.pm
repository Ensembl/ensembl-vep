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

use base qw(Bio::EnsEMBL::VEP::OutputFactory);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Utils qw(convert_arrayref);

my @OUTPUT_COLS = qw(
  Uploaded_variation
  Location
  Allele
  Gene
  Feature
  Feature_type
  Consequence
  cDNA_position
  CDS_position
  Protein_position
  Amino_acids
  Codons
  Existing_variation
);

my %OUTPUT_COLS_HASH = map {$_ => 1} @OUTPUT_COLS;

sub output_hash_to_line {
  my $self = shift;
  my $hash = shift;

  # "core" fields
  my @line = map {defined($hash->{$_}) ? convert_arrayref($hash->{$_}) : '-'} @OUTPUT_COLS;

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

sub field_order {
  my $self = shift;

  if(!exists($self->{field_order})) {

    my @extra_fields =
      map {@{$_->{fields}}}
      map {$_->[0]}
      grep {
        ref($_->[1]) eq 'ARRAY' ? scalar @{$_->[1]} : $_->[1]
      }
      map {[$_, $self->param($_->{flag})]}
      @{$self->flag_fields};
    
    $self->{field_order}->{$extra_fields[$_]} = $_ for 0..$#extra_fields;
  }

  return $self->{field_order};
}


return 1;
