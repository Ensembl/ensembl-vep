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

# EnsEMBL module for Bio::EnsEMBL::VEP::VariantRecoder
#
#

=head1 NAME

Bio::EnsEMBL::VEP::VariantRecoder - VariantRecoder runner class

=head1 SYNOPSIS

my $idt = Bio::EnsEMBL::VEP::VariantRecoder->new();
my $recoded = $idt->recode('rs699');

=head1 DESCRIPTION

The VariantRecoder class serves as a wrapper for a number of
VEP classes that is used to "recode" variant identifiers
to all possible alternatives:

- variant IDs
- HGVS genomic (g.)
- HGVS coding (c.)
- HGVS protein (p.)

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::VariantRecoder;

use base qw(Bio::EnsEMBL::VEP::Runner);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::VEP::Runner;
use Bio::EnsEMBL::VEP::Utils qw(find_in_ref merge_arrays);
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::Utils::VEP;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

=head2 new

  Arg 1      : hashref $config
  Example    : $runner = Bio::EnsEMBL::VEP::VariantRecoder->new($config);
  Description: Creates a new VariantRecoder object. The $config hash passed is
               used to create a Bio::EnsEMBL::VEP::Config object; see docs
               for this object and the variant_recoder script itself for allowed
               parameters.
  Returntype : Bio::EnsEMBL::VEP::VariantRecoder
  Exceptions : throws on invalid configuration, see Bio::EnsEMBL::VEP::Config
  Caller     : variant_recoder
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $config = shift || {};

  $config->{$_} = 1 for grep {!exists($config->{$_})} qw(
    database
    merged
    lrg
    check_existing
    failed
    no_prefetch
    hgvsg_use_accession
    ambiguous_hgvs
    no_stats
    json
    quiet
    buffer_size
  );

  $config->{fields} ||= 'id,hgvsg,hgvsc,hgvsp,spdi';

  # Flag to indicate it is a Variant Recoder job
  $config->{is_vr} = 1;

  my %set_fields = map {$_ => 1} ref($config->{fields}) eq 'ARRAY' ? @{$config->{fields}} : split(',', $config->{fields});

  # do some trickery to make sure we're not running unnecessary code
  # this first one only switches on the HGVS options for the requested fields  
  $config->{$_} = 1 for grep {$_ =~ /^hgvs/} keys %set_fields;
  $config->{$_} = 1 for grep {$_ =~ /^spdi/} keys %set_fields;

  # and this one switches on check_existing if the user wants variant IDs
  my %opt_map = ('id' => 'check_existing');
  $config->{$opt_map{$_}} = 1 for grep {$set_fields{$_}} keys %opt_map;

  # set up/down distance to 0, we only want overlaps
  $config->{distance} = 0;
  
  if($config->{vcf_string}){
    $config->{fields} = $config->{fields} . ',vcf_string';
  }

  if($config->{var_synonyms} && $config->{fields} !~ /var_synonyms/){
    $config->{fields} = $config->{fields} . ',var_synonyms';
  }

  my $self = $class->SUPER::new($config);

  return $self;
}


=head2 init

  Example    : $idt->init();
  Description: Runs some initialisation processes:
               - connect to DB
               - get annotation sources
               - internalise warnings
  Returntype : bool
  Caller     : recode(), recode_all()
  Status     : Stable

=cut

sub init {
  my $self = shift;

  return 1 if $self->{_initialized};

  $self->SUPER::init();

  $self->internalise_warnings();

  return 1;
}


=head2 recode_all

  Example    : my $results = $idt->recode_all();
  Description: Get all recoding results for the input file
               set up at initialisation or with $self->param('input_file')
  Returntype : hashref
  Caller     : variant_recoder
  Status     : Stable

=cut

sub recode_all {
  my $self = shift;

  $self->init();

  my $results = $self->_get_all_results();

  $self->finish();

  return $results;
}


=head2 recode
  
  Arg 1      : string $input_data
  Example    : my $results = $idt->recode('rs699');
  Description: Get recoding results for given input string
  Returntype : hashref
  Caller     : general
  Status     : Stable

=cut

sub recode {
  my $self = shift;
  my $input = shift;

  throw("ERROR: No input data supplied") unless $input;

  $self->param('input_data', $input);

  $self->init();
  my $results = $self->_get_all_results();

  $self->reset();

  return $results;
}


=head2 reset
  
  Example    : $idt->reset();
  Description: Reset input parameters. Used by recode() to prevent
               persistence of input, format etc settings between calls
               to recode()
  Caller     : recode()
  Status     : Stable

=cut

sub reset {
  my $self = shift;

  delete($self->{$_}) for qw(parser input_buffer);
  $self->param('format', 'guess');
  $self->param('input_data', undef);
}


=head2 _get_all_results
  
  Example    : my $results = $idt->_get_all_results();
  Description: Internal method used to fetch results for set up object.
  Caller     : recode(), recode_all()
  Status     : Stable

=cut

sub _get_all_results {
  my $self = shift;

  my $results = {};
  my $order   = [];

  my %want_keys = map {$_ => 1} @{$self->param('fields')};

  # Store data that is originally not linked to an allele
  my %keys_no_allele;
 
  if($want_keys{'id'}) {
    $keys_no_allele{'id'} = 1;
    delete($want_keys{'id'});
  }
  if($want_keys{'vcf_string'}) {
    $keys_no_allele{'vcf_string'} = 1;
    delete($want_keys{'vcf_string'});
  }
  if($want_keys{'var_synonyms'}) {
    $keys_no_allele{'var_synonyms'} = 1;
    delete($want_keys{'var_synonyms'});
  }

  while(my $line = $self->next_output_line(1)) {
    delete($line->{id});
    my $line_id = $line->{input};
   
    my %line_by_allele;
    my $consequences = $line->{transcript_consequences} ||= $line->{intergenic_consequences};

    # Split the consequences by alleles
    my %allele_consequence;
    foreach my $consequence (@$consequences) {
      my $allele = $consequence->{'variant_allele'};
      push @{$allele_consequence{$allele}}, $consequence;
    }

    # Hash contains all consequences separated by allele
    $line_by_allele{'consequences'} = \%allele_consequence;

    ##################
    ### vcf_string ###
    # Parse vcf string and build a hash by allele
    my %vcf_string_by_allele;
    if($keys_no_allele{'vcf_string'}) {
      # Minimise alleles - the alleles from vcf_string need to be compared with the alleles from allele_string
      # but we need to minimise both, otherwise we are comparing different things
      # example: rs1160446270 -> allele string is CCCCCCC but in vcf_string allele is CCC "22-16251519-CC-CCC"
      my $minimise_alleles_hash = $self->_minimise_allele($line->{'allele_string'}, $line->{start}, $line->{end}, $line->{strand}, 1);
      # Converts the allele_string to vcf format to compare with the original vcf_string, example: rs774003114
      my $allele_string_vcf_hash = $self->_convert_allele_to_vcf($line->{'allele_string'}, $line->{start}, $line->{end}, $line->{strand});

      # If there is more than one vcf_string then it's an array
      # example: rs56116432 -> "9-133256042-C-T" and "CHR_HG2030_PATCH-133256189-C-T"
      if(ref($line->{'vcf_string'})) {
        foreach my $vcf_string (@{$line->{'vcf_string'}}) {
            my @split_vcf = split /\-/, $vcf_string;
            my $alt_allele_from_vcf = $allele_string_vcf_hash->{$split_vcf[2].'/'.$split_vcf[3]};

            # Need to minimise the alleles for deletion and insertion
            if(!$alt_allele_from_vcf) {
              my $alt_allele_from_vcf_1 = $self->_minimise_allele_vcf($vcf_string, $line->{start}, $line->{end}, $line->{strand}, $minimise_alleles_hash, $allele_string_vcf_hash, 0);
              # If allele_string or vcf_string have reversed alleles then need to check the reverse for the comparison
              # Example: allele_string is 'G/C' and vcf_string is '21-25891856-C-G', 'G/C' can't be compared with 'C/G
              my $alt_allele_from_vcf_2 = $self->_minimise_allele_vcf($vcf_string, $line->{start}, $line->{end}, $line->{strand}, $minimise_alleles_hash, $allele_string_vcf_hash, 1);
            $alt_allele_from_vcf = $alt_allele_from_vcf_1 ? $alt_allele_from_vcf_1 : $alt_allele_from_vcf_2;
            }

            $vcf_string_by_allele{$alt_allele_from_vcf}->{'vcf_string'} = $vcf_string;
        }
      }
      else {
        # It only contains one vcf_string
          my @split_vcf = split /\-/, $line->{'vcf_string'};
          my $new_vcf_string = $split_vcf[2].'/'.$split_vcf[3];
          my $alt_allele_from_vcf = $allele_string_vcf_hash->{$split_vcf[2].'/'.$split_vcf[3]};

          if(!$alt_allele_from_vcf) {
            # minimise the alleles and check if it already exists (e.g. if it's the same as one of the alleles from allele_string or from allele_string converted to VCF)
            my $alt_allele_from_vcf_1 = $self->_minimise_allele_vcf($line->{'vcf_string'}, $line->{start}, $line->{end}, $line->{strand}, $minimise_alleles_hash, $allele_string_vcf_hash, 0);
            # Same as above but for the reverse alleles
            my $alt_allele_from_vcf_2 = $self->_minimise_allele_vcf($line->{'vcf_string'}, $line->{start}, $line->{end}, $line->{strand}, $minimise_alleles_hash, $allele_string_vcf_hash, 1);
            $alt_allele_from_vcf = $alt_allele_from_vcf_1 ? $alt_allele_from_vcf_1 : $alt_allele_from_vcf_2;
          }

          $vcf_string_by_allele{$alt_allele_from_vcf}->{'vcf_string'} = $line->{'vcf_string'};
      }
    }
    ### vcf_string ###
    ##################

    ##################
    ####### ID #######
    # Parse ID and build a hash by allele
    # COSMIC and HGMD IDs don't have an allele - they are annexed to all alleles
    if($keys_no_allele{'id'}) {

      my @ids_no_allele; # stores the cosmic and hgms ids
      foreach my $co_var (@{$line->{'colocated_variants'}}) {
        # Store the COSMIC and HGMD IDs (not linked to an allele) for later
        if($co_var->{'allele_string'} =~ /COSMIC|HGMD/) {
          push @ids_no_allele, $co_var->{'id'}; 
        }
        else {

          my @split_allele = split /\//, $co_var->{'allele_string'};
          # delete ref allele - we only want to check the alt alleles from the colocated variants
          shift @split_allele;
          foreach my $allele (@split_allele) {
            if($allele_consequence{$allele}) {
              push @{$vcf_string_by_allele{$allele}->{'id'}}, $co_var->{'id'};
            }
            else {
              # If the allele is not stored then we need to check for the reverse
              my $allele_rev = $allele;
              reverse_comp(\$allele_rev);
              # We also check the minimised allele
              my $minimise_co_allele = $self->_minimise_allele($co_var->{'allele_string'}, $line->{start}, $line->{end}, $line->{strand}, 0);
              $minimise_co_allele =~ s/.*\///;
              my $xxx = $allele_consequence{$allele_rev} ? $allele_rev : $minimise_co_allele;

              push @{$vcf_string_by_allele{$xxx}->{'id'}}, $co_var->{'id'};
            }
          }
        }
      }
      # Link the COSMIC and HGMD IDs to all alleles
      if(scalar(@ids_no_allele) != 0) { 
        foreach my $key_allele (keys %{$line_by_allele{'consequences'}}) {
          foreach my $id_no_allele (@ids_no_allele) {
            push @{$vcf_string_by_allele{$key_allele}->{'id'}}, $id_no_allele;
          }
        }
      }
    }
    ####### ID #######
    ##################

    ################################
    ####### Variant synonyms #######
    # Attach variant synonyms to hash by allele
    if($line->{'var_synonyms'} && $keys_no_allele{'var_synonyms'}) {
      foreach my $key_allele (keys %{$line_by_allele{'consequences'}}) {
      $vcf_string_by_allele{$key_allele}->{'var_synonyms'} = $line->{'var_synonyms'};
      }
    }

    merge_arrays($order, [$line_id]);

    foreach my $allele (keys %{$line_by_allele{'consequences'}}) {
      find_in_ref($line_by_allele{'consequences'}->{$allele}, \%want_keys, $results->{$line_id}->{$allele} ||= {input => $line_id});
      find_in_ref($vcf_string_by_allele{$allele}, \%keys_no_allele, $results->{$line_id}->{$allele} ||= {input => $line_id});
    }

    if(@{$self->warnings}) {
      $results->{$line_id}->{warnings} = [map {$_->{msg}} @{$self->warnings}];
      $self->internalise_warnings();
    }
  }

  return [map {$results->{$_}} @$order];
}

=head2 _minimise_allele

  Example    : my $allele = $idt->_minimise_allele($allele_string, $start, $end, $strand);
  Description: Internal method used to minimise the alleles for comparison.
  Caller     : _get_all_results()
  Status     : Stable

=cut

sub _minimise_allele {
  my ($self, $allele_string, $start, $end, $strand, $flag) = @_;

  my $minimise_alleles;
  my $result_allele;

  my @split_allele = split /\//, $allele_string;
  my $ref_allele = shift @split_allele;

  foreach my $alt_allele (@split_allele) {
    my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
        -start   => $start,
        -end     => $end,
        -strand  => $strand,
        -allele_string => $ref_allele.'/'.$alt_allele,
    );

    my $minimise_vf = Bio::EnsEMBL::Variation::Utils::VEP->minimise_alleles([$vf]);
    my $final_allele = $minimise_vf->[0]->{allele_string};

    $minimise_alleles->{$final_allele} = $alt_allele;

    if(!$flag) {
      $result_allele = $final_allele;
    }
  }

  return $flag ? $minimise_alleles : $result_allele;
}

=head2 _minimise_allele_vcf

  Example    : my $allele = $idt->_minimise_allele_vcf($vcf_string, $start, $end, $strand, $hash);
  Description: Internal method used to minimise the VCF alleles for comparison.
  Caller     : _get_all_results()
  Status     : Stable

=cut

sub _minimise_allele_vcf {
  my ($self, $vcf_string, $start, $end, $strand, $minimise_alleles_hash, $allele_string_vcf_hash, $reverse) = @_;

  my @split_vcf = split /\-/, $vcf_string;
  # Example: allele_string is 'G/C' and vcf_string is '21-25891856-C-G', 'G/C' can't be compared with 'C/G'
  # If reverse = 1 then the new vcf_string is 'G/C' which can be compared with allele_string 'G/C'
  my $ref_allele_vcf = $reverse == 0 ? $split_vcf[2] : $split_vcf[-1];
  my $alt_allele_vcf = $reverse == 0 ? $split_vcf[-1] : $split_vcf[2];

  my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
      -start   => $start,
      -end     => $end,
      -strand  => $strand,
      -allele_string => $ref_allele_vcf.'/'.$alt_allele_vcf,
  );

  my $minimise_vf = Bio::EnsEMBL::Variation::Utils::VEP->minimise_alleles([$vf]);
  my $after_minimise_alleles = $minimise_vf->[0]->{allele_string};
  my $final_alt_allele = defined ($minimise_alleles_hash->{$after_minimise_alleles}) ? $minimise_alleles_hash->{$after_minimise_alleles} : $allele_string_vcf_hash->{$after_minimise_alleles};

  return $final_alt_allele;
}

=head2 _convert_allele_to_vcf

  Example    : my $hash = $idt->_convert_allele_to_vcf($allele_string, $start, $end, $strand);
  Description: Internal method used to convert an allele string to the VCF format.
               Returns a hash mapping the VCF format with the allele string.
  Caller     : _get_all_results()
  Status     : Stable

=cut

sub _convert_allele_to_vcf {
  my ($self, $allele_string, $start, $end, $strand) = @_;

  my $vcf_alleles;

  my @split_allele = split /\//, $allele_string;
  my $ref_allele = shift @split_allele;

  foreach my $alt_allele (@split_allele) {
    my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
        -start   => $start,
        -end     => $end,
        -strand  => $strand,
        -allele_string => $ref_allele.'/'.$alt_allele,
    );

    my $converted_to_vcf = $vf->to_VCF_record;
    my $ref_allele_vcf = ${$converted_to_vcf}[3];
    my $alt_allele_vcf = ${$converted_to_vcf}[4];

    $ref_allele_vcf =~ s/N//;
    $alt_allele_vcf =~ s/N//;
    if(!$ref_allele_vcf) {
      $ref_allele_vcf = '-';
    }
    if(!$alt_allele_vcf) {
      $alt_allele_vcf = '-';
    }

    my @split_vcf_allele = split /,/, $alt_allele_vcf;
    foreach my $alt_split (@split_vcf_allele) {
      $vcf_alleles->{$ref_allele_vcf.'/'.$alt_split} = $alt_allele;
    }
  }

  return $vcf_alleles;
}

1;
