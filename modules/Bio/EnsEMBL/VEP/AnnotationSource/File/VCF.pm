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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::VCF - VCF annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::VCF;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::VEP::Utils qw(trim_sequences);
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  $self->fields($hashref->{fields}) if $hashref->{fields};

  return $self;
}

sub fields {
  my $self = shift;
  $self->{fields} = shift if @_;
  return $self->{fields};
}

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($self->file);
}

sub _create_records {
  my $self = shift;
  my $overlap_result = shift;

  my @records;

  # get the fields data from the parser
  my $fields_data = {};

  if(my $fields = $self->fields) {
    my $parser = $self->parser;
    my $info = $parser->get_info;

    foreach my $field(grep {exists($info->{$_})} @$fields) {
      my $value = $info->{$field};
      my @return;

      if($value =~ /\,/) {
        my @split = split(',', $value);
        $fields_data->{$field}->{$_} = $split[$_] for 0..$#split;
      }
      else {
        $fields_data->{$field} = $value;
      }
    }
  }

  # exact match returns a hashref
  if(ref($overlap_result) eq 'HASH') {

    while(my $allele = shift @{$overlap_result->{alleles}}) {
      my $index = shift @{$overlap_result->{indexes}};

      my $record = {
        name => $self->_get_record_name,
        allele => $allele
      };

      foreach my $field(keys %$fields_data) {
        my $data = $fields_data->{$field};
        $record->{fields}->{$field} = ref($data) eq 'HASH' ? $data->{$index} : $data;
      }

      push @records, $record;
    }
  }
  else {
    my $record = {
      name => $self->_get_record_name,
    };

    foreach my $field(keys %$fields_data) {
      my $data = $fields_data->{$field};

      if(ref($data) eq 'HASH') {
        $record->{fields}->{$field} = join(',', map {$data->{$_}} sort {$a <=> $b} keys %{$fields_data->{$field}});
      }
      else {
        $record->{fields}->{$field} = $data;
      }
    }

    push @records, $record;
  }

  return \@records;
}

sub _get_record_name {
  my $self = shift;
  my $parser = $self->parser;

  my $id = $parser->get_IDs->[0];

  return ($self->report_coords || !$id || $id eq '.') ?
    sprintf(
      '%s:%i-%i',
      $parser->get_seqname,
      $parser->get_start,
      $parser->get_end
    ) :
    $id;
}

sub _record_overlaps_VF {
  my $self = shift;
  my ($vf) = @_;

  # we can use the superclass method if overlap type
  return $self->SUPER::_record_overlaps_VF(@_)
    if $self->type eq 'overlap' || ref($vf) eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature';

  # exact more difficult, we need to check each allele
  my $parser = $self->parser;

  my $vf_ref  = $vf->ref_allele_string;
  my %vf_alts = map {$_ => 1} @{$vf->alt_alleles};
  my $vf_strand = $vf->strand;

  # use these as backups as we might modify them
  my $orig_start = $parser->get_raw_start;
  my $orig_ref   = $parser->get_reference;

  # we're going to return a hashref containing the matched allele indexes and alleles
  # that way it also passes a boolean check
  my @matched_indexes;
  my @matched_vf_alleles;
  my $i = 0;
  foreach my $alt(@{$parser->get_alternatives}) {

    # we're going to minimise the VCF alleles one by one and compare the coords
    # first make copies of everything
    my $start = $orig_start;
    my $ref   = $orig_ref;

    ($ref, $alt, $start) = @{trim_sequences($ref, $alt, $start)};
    $ref ||= '-';
    $alt ||= '-';

    # check strand before doing expensive revcomp
    if($vf->{start} == $start) {
      if($vf_strand < 0) {
        reverse_comp(\$ref);
        reverse_comp(\$alt);
      }

      if($vf_ref eq $ref && $vf_alts{$alt}) {
        push @matched_indexes, $i;
        push @matched_vf_alleles, $alt;
      }
    }
    $i++;
  }

  return @matched_indexes ? {indexes => \@matched_indexes, alleles => \@matched_vf_alleles} : 0;
}


1;