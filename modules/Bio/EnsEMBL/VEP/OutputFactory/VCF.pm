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

# EnsEMBL module for Bio::EnsEMBL::VEP::OutputFactory::VCF
#
#

=head1 NAME

Bio::EnsEMBL::VEP::OutputFactory::VCF - VCF format output factory

=head1 SYNOPSIS

my $of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({
  config => $config,
});

# print headers
print "$_\n" for @{$of->headers};

# print output
print "$_\n" for @{$of->get_all_lines_by_InputBuffer($ib)};

=head1 DESCRIPTION

An OutputFactory class to generate VCF-format output.

If the input was VCF, the original input line is appended with
VEP annotations.

If the input was another format, VCF lines are generated from
scratch. In the case of unbalanced substitutions (e.g. insertions
or deletions), this means that the base preceding the variant
must be prepended to the REF and ALT alleles to comply with 
VCF specification.

The VEP output data itself is concatenated in a carefully delimited
string in the INFO field of the VCF, using the (configurable) key
CSQ. The value for the CSQ field consists of one or more
comma-separated "chunks" of annotation representing
one variant allele/feature combination. Each chunk contains a
constant number of fields, separated by pipe ("|") characters, with
the fields being defined in the VCF header. Empty values for fields
are represented by empty strings, meaning multiple pipe characters
may be found consecutively.

The nested nature of this format means that several character
substitutions must be used within VEP fields to preserve VCF format:

,          ==> &
=          ==> %3B
whitespace ==> _

=head1 METHODS


=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::OutputFactory::VCF;

use base qw(Bio::EnsEMBL::VEP::OutputFactory);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::VEP::Utils qw(convert_arrayref get_version_data);
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


=head2 new

  Arg 1      : hashref $args
  Example    : $of = Bio::EnsEMBL::VEP::OutputFactory::VCF->new({
                 config => $config,
               });
  Description: Creates a new Bio::EnsEMBL::VEP::OutputFactory::VCF object.
               Has its own constructor to add several params via
               add_shortcuts()
  Returntype : Bio::EnsEMBL::VEP::OutputFactory::VCF
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(
    fields
    vcf_info_field
    keep_csq
    web_output
  )]);

  $self->{input_format} = $self->param('format');

  return $self;
}


=head2 headers

  Example    : $headers = $of->headers();
  Description: Get list of headers to print out.
  Returntype : listref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub headers {
  my $self = shift;

  my $info = $self->header_info;
  my $field_descs = \%Bio::EnsEMBL::VEP::Constants::FIELD_DESCRIPTIONS;

  # VCFs have metadata headers starting with ##
  # and one line of column headers starting with #
  my (@headers, $col_heading);
  
  my $vcf_fileformat = '##fileformat=VCFv4.1';

  # input was VCF
  if($info->{input_headers} && scalar @{$info->{input_headers}}) {

    # If the VCF file format header was missing in the input VCF file
    if (!grep { $_ =~ /^##fileformat=VCFv4/i } @{$info->{input_headers}}) {
      @headers = ($vcf_fileformat);
    }

    if ($self->{keep_csq}) {
      push @headers, @{$info->{input_headers}};
    }
    else {
      my $fieldname = $self->{vcf_info_field} || 'CSQ';
      foreach my $input_header (@{$info->{input_headers}}) {
        push @headers, $input_header unless ($input_header =~ /ID=$fieldname,/i || $input_header =~ /^##VEP/);
      }
    }
    $col_heading = pop @headers;
  }

  # input wasn't VCF
  else {
    @headers = ($vcf_fileformat);
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

  # add API module info
  my $software_version_data = get_version_data();
  my $software_version_string = join(' ',
    map {
      sprintf(
        '%s=%s%s',
        $_,
        $software_version_data->{$_}->{release},
        (defined($software_version_data->{$_}->{sub}) ? '.'.substr($software_version_data->{$_}->{sub}, 0, 7) : '')
      )
    }
    grep {$_ ne 'ensembl-vep'} keys %{$software_version_data}
  );
  $headers[-1] .= ' '.$software_version_string if $software_version_string;

  # add misc version data
  $headers[-1] .= ' '.join(' ',
    map {$_.'="'.$info->{version_data}->{$_}.'"'}
    grep {defined($info->{version_data}->{$_})}
    sort keys %{$info->{version_data}}
  ) if $info->{version_data};

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


=head2 get_all_lines_by_InputBuffer

  Arg 1      : Bio::EnsEMBL::VEP::InputBuffer $ib
  Example    : $lines = $of->get_all_lines_by_InputBuffer($ib);
  Description: Gets all lines (strings suitable for writing to output) given
               an annotated input buffer, one line per input variant.
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub get_all_lines_by_InputBuffer {
  my $self = shift;
  my $buffer = shift;

  my @return;

  map {@{$self->reset_shifted_positions($_)}}
    @{$buffer->buffer};

  $self->rejoin_variants_in_InputBuffer($buffer) if $buffer->rejoin_required;

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
      if($line->[7] =~ /(^|\;)$fieldname\=/ && !$self->{keep_csq}) {
        $line->[7] =~ s/(^|\;)$fieldname\=\S+?(\;\S|$)/$2/;
        $line->[7] =~ s/^\;//;
      }
    }

    # we have to create one if input wasnt VCF
    else {
      $vf->{slice} ||= $self->get_slice($vf->{chr});
      $line = $vf->to_VCF_record();
      $line->[7] = '' if ($line->[7] || '') eq '.';
    }

    my @chunks =
      map {$self->output_hash_to_vcf_info_chunk($_, $vf->strand)} 
      @{$self->get_all_output_hashes_by_VariationFeature($vf)};

    if(@chunks) {
      $line->[7] .= ';' if $line->[7];
      $line->[7] .= $fieldname.'='.join(",", @chunks);
    }
    else {
      $line->[7] ||= '.';
    }

    push @return, join("\t", map {defined($_) ? $_ : '.'} @$line);

    $self->write_web_output($vf) if $self->{web_output};
  }

  return \@return;
}


=head2 output_hash_to_vcf_info_chunk

  Arg 1      : hashref $vf_hash
  Arg 2      : (optional) int $vf_strand
  Example    : $chunk = $of->output_hash_to_vcf_info_chunk($vf_hash);
  Description: Converts a hashref as retrieved from
               get_all_output_hashes_by_VariationFeature to a pipe-separated
               chunk suitable for adding to the CSQ key in the VCF INFO field.
  Returntype : string
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

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

      if(defined($data)) {
        $data =~ s/\,/\&/g;
        $data =~ s/\;/\%3B/g;
        $data =~ s/\s+/\_/g;
        $data =~ s/\|/\&/g;
      }

      push @chunk, $data;
    }
    # keep 0 values
    elsif (defined $hash->{$col} && $hash->{$col} ne '' 
           && ($hash->{$col} == 0)) {
      push @chunk, $hash->{$col};
    }
    else {
      push @chunk, '';
    }
  }

  return join('|', @chunk);
}


=head2 fields

  Example    : $fields = $of->fields();
  Description: Gets list of fields to be populated
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : output_hash_to_vcf_info_chunk()
  Status     : Stable

=cut

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


=head2 write_web_output

  Arg 1      : Bio::EnsEMBL::Variation::BaseVariationFeature $bvf
  Example    : $of->write_web_output($bvf);
  Description: Writes a line of summarised output to the web_output
               filehandle.

               The VEP web backend uses VCF as its output, but restrictions
               in filesystems means that VEP must write a secondary
               output file containing summary data to allow user input
               to be rendered as a track on region-in-detail view etc.
  Returntype : string
  Exceptions : none
  Caller     : get_all_lines_by_InputBuffer(), 
  Status     : Stable

=cut

sub write_web_output {
  my $self = shift;
  my $vf = shift;

  my $fh = $self->web_output_fh;

  my $as = $vf->{allele_string};
  my $id = $vf->{variation_name};
  if(defined($as) && length($as) > 50) {
    my @new_alleles;
    
    foreach my $allele(split(/\//, $as)) {
      if(length($allele) > 50) {
        my $new = length($allele).'BP_SEQ';
        push @new_alleles, $new;
        
        $id =~ s/$allele/$new/e;
      }
      else {
        push @new_alleles, $allele;
      }
    }
    
    $as = join("/", @new_alleles);
  }

  printf $fh "%s\t%i\t%i\t%s\t%s\t%s\t%s\n",
    $vf->{chr}, $vf->{start}, $vf->{end},
    $as || $vf->class_SO_term, 1,
    $id,
    $vf->display_consequence;
}


=head2 web_output_fh

  Example    : $fh = $of->web_output_fh();
  Description: Gets filehandle for writing web output to.
  Returntype : glob
  Exceptions : throws if unable to write to file
  Caller     : write_web_output(), 
  Status     : Stable

=cut

sub web_output_fh {
  my $self = shift;

  if(!exists($self->{_web_output_fh})) {
    my $fh = FileHandle->new();
    $fh->open(">".$self->{web_output}) or throw $!;
    $self->{_web_output_fh} = $fh;
  }

  return $self->{_web_output_fh};
}

1;
