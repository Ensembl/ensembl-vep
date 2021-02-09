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

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::VCF - VCF annotation source

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::VCF->new({
  config => $config,
  file   => "my_variants.vcf.gz",
  type   => "exact"
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

VCF format custom annotation source. VCF files must be chromosome/pos
sorted, compressed with bgzip and indexed with tabix.

Additional fields of data from the VCF INFO field may be added to the
custom annotation data returned by setting fields().

Using the "exact" annotation type with this class enforces allele as
well as positional checking of the data, and any INFO fields requested
that have allele-specific data will have only the relevant data added
depending on the alleles of the variant from the user-fed InputBuffer.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::VCF;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File);

=head2 new

  Arg 1      : hashref $args
               {
                 config        => Bio::EnsEMBL::VEP::Config $config,
                 file          => string $filename,
                 short_name    => (optional) string $short_name,
                 type          => (optional) string $type (overlap (default), exact),
                 report_coords => (optional) bool $report_coords,
                 fields        => arrayref $INFO_fields_to_add
               }
  Example    : $as = Bio::EnsEMBL::VEP::AnnotationSource::File::VCF->new($args);
  Description: Create a new Bio::EnsEMBL::VEP::AnnotationSource::File::VCF object.
  Returntype : Bio::EnsEMBL::VEP::AnnotationSource::File::VCF
  Exceptions : none
  Caller     : Bio::EnsEMBL::VEP::AnnotationSource::File->new()
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  my $hashref = $_[0];

  $self->add_shortcuts(['custom_multi_allelic']);

  $self->fields($hashref->{fields}) if $hashref->{fields};      ## report INFO & FILTER fields

  return $self;
}


=head2 fields

  Arg 1      : (optional) arrayref $fields
  Example    : $fields = $as->fields()
  Description: Getter/setter for the list of fields to be added from the VCF INFO field.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fields {
  my $self = shift;
  $self->{fields} = shift if @_;
  return $self->{fields};
}

=head2 parser

  Example    : $parser = $as->parser();
  Description: Get ensembl-io parser to read from file
  Returntype : Bio::EnsEMBL::IO::Parser::VCF4Tabix
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub parser {
  my $self = shift;
  return $self->{parser} ||= Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($self->file);
}


=head2 _create_records
 
  Arg 1      : bool or hashref $overlap_result
  Example    : $records = $as->_create_records($overlap_result);
  Description: Create a custom annotation record from the current
               record as read from the annotation source. Adds INFO field
               data if requested by setting fields().
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : annotate_VariationFeature()
  Status     : Stable

=cut

sub _create_records {
  my $self = shift;
  my $overlap_result = shift;

  my @records;

  # get the fields data from the parser
  my $fields_data = {};

  if(my $fields = $self->fields) {
    my $parser = $self->parser;
    my $info = $parser->get_info;
    my $alt_count = scalar @{$parser->get_alternatives};

    foreach my $field(grep {exists($info->{$_})} @$fields) {
      my $value = $info->{$field};

      # substitute in a value of 1 for absent; some keys in VCF do not have a value
      # but VEP needs something to write out
      $value = 1 if !defined($value);
      
      my @return;
      my $metadata = $parser->{metadata} || {};
      my $is_source_clinvar = defined($metadata->{source}) ? $metadata->{source} eq 'ClinVar' : 0;	
      
      if(!($self->{'custom_multi_allelic'} || $is_source_clinvar) && $value =~ /\,/) {
        my @split = split(',', $value);

        # some VCFs have data for REF included
        shift @split if scalar(@split) == $alt_count + 1;

        # don't do allele-specific if @split now doesn't match alt_count
        if(scalar(@split) == $alt_count) {
          $fields_data->{$field}->{$_} = $split[$_] for 0..$#split;  
        }
        else {
          $fields_data->{$field} = $value;  
        }
        
      }
      else {
        $fields_data->{$field} = $value;
      }
    }
    ## extract pass/fail info from filter column
    $fields_data->{FILTER} .= $parser->get_raw_filter_results();
    $fields_data->{FILTER} =~ s/\;/\,/g;
  }


  # exact match returns a arrayref
  if(ref($overlap_result) eq 'ARRAY') {

    foreach my $result_hash(@$overlap_result) {
      my $index  = $result_hash->{b_index};

      my $record = {
        name => $self->_get_record_name,
        allele => $result_hash->{a_allele}
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


=head2 _get_record_name
 
  Example    : $record_name = $as->_get_record_name();
  Description: Get name for the current record using either ID as
               found in the source or record coordinates. Defaults
               to using first ID if multiple found, or coordinates
               if none found.
  Returntype : string
  Exceptions : none
  Caller     : _create_records()
  Status     : Stable

=cut

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


=head2 _record_overlaps_VF
 
  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature
  Example    : $overlaps = $as->_record_overlaps_VF($vf);
  Description: Determine whether the given VariationFeature overlaps
               the current record, depending on the set type(). Note
               for "exact" type allele as well as position is matched.
               VCF entries with multiple ALTs are treated as separate
               REF/ALT pairs, with the REF/ALT being trimmed to the minimum
               shared sequence before comparison to the alleles from
               the user input variant.

               Returns an arrayref of hashrefs representing the matching
               alleles and their indexes:
               [
                 {
                   a_index  => $VF_index_1,
                   a_allele => $VF_allele_1,
                   b_index  => $VCF_ALT_index_1,
                   b_allele => $VCF_allele_1,
                 },
                 {
                   a_index  => $VF_index_2,
                   a_allele => $VF_allele_2,
                   b_index  => $VCF_ALT_index_2,
                   b_allele => $VCF_allele_2,
                 },
               ]
  Returntype : arrayref
  Exceptions : none
  Caller     : annotate_VariationFeature()
  Status     : Stable

=cut

sub _record_overlaps_VF {
  my $self = shift;
  my ($vf) = @_;

  # we can use the superclass method if overlap type
  return $self->SUPER::_record_overlaps_VF(@_)
    if $self->type eq 'overlap' || ref($vf) eq 'Bio::EnsEMBL::Variation::StructuralVariationFeature';

  # exact more difficult, we need to check each allele
  my $parser = $self->parser;

  my $matches = get_matched_variant_alleles(
    {
      ref    => $vf->ref_allele_string,
      alts   => $vf->alt_alleles,
      pos    => $vf->{start},
      strand => $vf->strand
    },
    {
      ref  => $parser->get_reference,
      alts => $parser->get_alternatives,
      pos  => $parser->get_raw_start,
    }
  );

  return @$matches ? $matches : 0;
}


1;
