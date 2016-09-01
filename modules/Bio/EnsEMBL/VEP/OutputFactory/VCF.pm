=head1 LICENSE

Copyright [2016] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::OutputFactory::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::OutputFactory::VCF - VCF format output factory

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::OutputFactory::VCF;

use base qw(Bio::EnsEMBL::VEP::OutputFactory);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::VEP::Utils qw(convert_arrayref);
use Bio::EnsEMBL::VEP::Constants;

my @VCF_COLS = qw(
  Allele
  Consequence
  IMPACT
  SYMBOL
  Gene
  Feature_type
  Feature
  BIOTYPE
  EXON
  INTRON
  HGVSc
  HGVSp
  cDNA_position
  CDS_position
  Protein_position
  Amino_acids
  Codons
  Existing_variation
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    fields
    vcf_info_field
    keep_csq
  )]);

  $self->{input_format} = $self->param('format');

  return $self;
}

sub headers {
  my $self = shift;

  my $info = $self->header_info;
  my $field_descs = \%Bio::EnsEMBL::VEP::Constants::FIELD_DESCRIPTIONS;

  # VCFs have metadata headers starting with ##
  # and one line of column headers starting with #
  my (@headers, $col_heading);
  
  # input was VCF
  if($info->{input_headers} && scalar @{$info->{input_headers}}) {
    my @input_headers = @{$info->{input_headers}};

    my $col_heading_pair = pop @input_headers;
    $col_heading = '#'.join("\t", @{$col_heading_pair->[1]});

    # add the remaining headers
    @headers = map {sprintf('##%s=%s', $_->[0], $_->[1])} @input_headers;
  }

  # input wasn't VCF
  else {
    @headers = ('##fileformat=VCFv4.1');
    $col_heading = '#'.join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO));
  }

  # add VEP version string
  push @headers, sprintf(
    '##VEP="v%i" time="%s"%s%s',
    $info->{vep_version},
    $info->{time},
    $info->{cache_dir} ? ' cache="'.$info->{cache_dir}.'"' : '',
    $info->{db_name} ? ' db="'.$info->{db_name}.'@'.$info->{db_host}.'"' : ''
  );

  # add misc version data
  $headers[-1] .= ' '.join(' ', map {$_.'="'.$info->{version_data}->{$_}.'"'} sort keys %{$info->{version_data}}) if $info->{version_data};

  # add VEP column def
  push @headers, sprintf(
    '##INFO=<ID=%s,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: %s">',
    $self->{vcf_info_field},
    join("|", @{$self->fields})
  );

  # plugin headers
  push @headers, map {sprintf('##%s=%s', $_->[0], $_->[1])} @{$self->get_plugin_headers};

  # custom headers
  push @headers, map {sprintf('##INFO=<ID=%s,Number=.,Type=String,Description="%s">', $_->[0], $_->[1])} @{$self->get_custom_headers};

  push @headers, $col_heading;

  return \@headers;
}

sub get_all_lines_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my @return;

  foreach my $vf(@{$buffer->buffer}) {

    my $line;
    my $fieldname = $self->{vcf_info_field} || 'CSQ';

    # if input was VCF then we get _line with the original contents
    if($self->{input_format} eq 'vcf' && $vf->{_line}) {

      $line = $vf->{_line};

      # array copy to keep original intact
      my @tmp = @$line;
      $line = \@tmp;

      if(!defined($line->[7]) || $line->[7] eq '.') {
        $line->[7] = '';
      }

      # nuke existing CSQ field?
      if($line->[7] =~ /$fieldname\=/ && !$self->{keep_csq}) {
        $line->[7] =~ s/$fieldname\=\S+?(\;|$)(\S|$)/$2/;
      }
    }

    # we have to create one if input wasnt VCF
    else {
      $line = $self->VariationFeature_to_VCF_record($vf);
      $line->[7] = '';
    }

    $line->[7] .= ';' if $line->[7];

    $line->[7] .= $fieldname.'='.join(",",
      map {$self->output_hash_to_vcf_info_chunk($_, $vf->strand)} 
      @{$self->get_all_output_hashes_by_VariationFeature($vf)}
    );

    push @return, join("\t", @$line);
  }

  return \@return;
}


sub output_hash_to_vcf_info_chunk {
  my $self = shift;
  my $hash = shift;
  my $strand = shift || 1;

  my @chunk;

  # use the field list (can be user-defined by setting --fields)
  for my $col(@{$self->fields}) {

    # search for data in main line hash as well as extra field
    if(my $data = $hash->{$col}) {
      $data = convert_arrayref($data, '&');

      if($col eq 'Allele') {
        reverse_comp(\$data) if $strand < 0;
      }
      else {
        $data = '' if $data eq '-';  
      }
      
      $data =~ s/\;/\%3B/g if defined $data;

      push @chunk, $data;
    }
    else {
      push @chunk, '';
    }
  }

  return join('|', @chunk);
}

sub fields {
  my $self = shift;

  if(!defined($self->{fields})) {

    my @fields = @VCF_COLS;

    my %vcf_cols = map {$_ => 1} @VCF_COLS;
    
    @fields = @VCF_COLS;
    
    push @fields, 
      grep {!$vcf_cols{$_}}
      @{$self->flag_fields};

    push @fields, map {$_->[0]} @{$self->get_plugin_headers}, @{$self->get_custom_headers};
    
    $self->{fields} = \@fields;
  }

  return $self->{fields};
}

sub VariationFeature_to_VCF_record {
  my $self = shift;
  my $vf = shift;

  # look for imbalance in the allele string
  if(ref($vf) eq 'Bio::EnsEMBL::Variation::VariationFeature') {
    my %allele_lengths;
    my @alleles = split '\/', $vf->allele_string;

    map {reverse_comp(\$_)} @alleles if $vf->strand < 0;

    foreach my $allele(@alleles) {
      $allele =~ s/\-//g;
      $allele_lengths{length($allele)} = 1;
    }

    # in/del/unbalanced
    if(scalar keys %allele_lengths > 1) {

      my $prev_base = $self->get_prev_base($vf);

      for my $i(0..$#alleles) {
        $alleles[$i] =~ s/\-//g;
        $alleles[$i] = $prev_base.$alleles[$i];
      }

      return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start - 1,
        $vf->variation_name || '.',
        shift @alleles,
        (join ",", @alleles),
        '.', '.', '.'
      ];

    }

    # balanced sub
    else {
      return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start,
        $vf->variation_name || '.',
        shift @alleles,
        (join ",", @alleles),
        '.', '.', '.'
      ];
    }
  }

  # SV
  else {

    # convert to SO term
    my %terms = (
      insertion => 'INS',
      deletion => 'DEL',
      tandem_duplication => 'TDUP',
      duplication => 'DUP'
    );

    my $alt = '<'.($terms{$vf->class_SO_term} || $vf->class_SO_term).'>';

    return [
      $vf->{chr} || $vf->seq_region_name,
      $vf->start - 1,
      $vf->variation_name || '.',
      $self->get_prev_base($vf),
      $alt,
      '.', '.',
      'END='.$vf->end
    ];
  }
}

sub get_prev_base {
  my $self = shift;
  my $vf = shift;  

  # we need the ref base before the variation
  # default to N in case we cant get it
  my $prev_base = 'N';

  $vf->{slice} ||= $self->get_slice($vf->{chr});

  if(defined($vf->{slice})) {
    my $slice = $vf->slice->sub_Slice($vf->start - 1, $vf->start - 1);
    $prev_base = $slice->seq if defined($slice);
  }

  return $prev_base;
}

1;