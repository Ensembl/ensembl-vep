=head1 LICENSE

Copyright [2016-2020] EMBL-European Bioinformatics Institute

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
use Data::Dumper;

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

    $line_by_allele{'consequences'} = \%allele_consequence;

    # Parse vcf string and build a hash by allele
    my %vcf_string_by_allele;
    if($keys_no_allele{'vcf_string'}) {
      # Minimise alleles
      my $minimise_alleles_hash = $self->_minimise_allele($line->{'allele_string'}, $line->{start}, $line->{end}, $line->{strand});
      # If there is more than one vcf_string then it's an array
      if(ref($line->{'vcf_string'})) {
        foreach my $vcf_string (@{$line->{'vcf_string'}}) {
          my $alt_allele_vcf_2 = $self->_minimise_allele_vcf($vcf_string, $line->{start}, $line->{end}, $line->{strand}, $minimise_alleles_hash);
          $vcf_string_by_allele{$alt_allele_vcf_2}->{'vcf_string'} = $vcf_string;
        }
      }
      else {
        my $alt_allele_vcf = $self->_minimise_allele_vcf($line->{'vcf_string'}, $line->{start}, $line->{end}, $line->{strand}, $minimise_alleles_hash);
        $vcf_string_by_allele{$alt_allele_vcf}->{'vcf_string'} = $line->{'vcf_string'};
      }
    }

    # Parse ID and build a hash by allele
    # COSMIC and HGMD IDs don't have an allele - they are annexed to all alleles
    if($keys_no_allele{'id'}) {
      my @ids_no_allele;
      foreach my $co_var (@{$line->{'colocated_variants'}}) {
        # Store the COSMIC and HGMD IDs (not linked to an allele) for later 
        if($co_var->{'allele_string'} =~ /COSMIC|HGMD/) {
          push @ids_no_allele, $co_var->{'id'}; 
        }
        else {
          my @split_allele = split /\//, $co_var->{'allele_string'};
          # delete ref allele
          shift @split_allele;
          foreach my $allele (@split_allele) {
            push @{$vcf_string_by_allele{$allele}->{'id'}}, $co_var->{'id'};
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
  Caller     : recode(), recode_all()
  Status     : Stable

=cut

sub _minimise_allele {
  my ($self, $allele_string, $start, $end, $strand) = @_;

  my $minimise_alleles;

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
  }

  return $minimise_alleles;
}

sub _minimise_allele_vcf {
  my ($self, $vcf_string, $start, $end, $strand, $minimise_alleles_hash) = @_;

  my @split_vcf = split /\-/, $vcf_string;
  my $ref_allele_vcf = $split_vcf[2];
  my $alt_allele_vcf = $split_vcf[-1];

  my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
      -start   => $start,
      -end     => $end,
      -strand  => $strand,
      -allele_string => $ref_allele_vcf.'/'.$alt_allele_vcf,
  );

  my $minimise_vf = Bio::EnsEMBL::Variation::Utils::VEP->minimise_alleles([$vf]);
  my $after_minimise_alleles = $minimise_vf->[0]->{allele_string};
  my $final_alt_allele = $minimise_alleles_hash->{$after_minimise_alleles};

  return $final_alt_allele;
}

1;
