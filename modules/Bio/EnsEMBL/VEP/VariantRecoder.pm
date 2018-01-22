=head1 LICENSE

Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

  $config->{fields} ||= 'id,hgvsg,hgvsc,hgvsp';

  my %set_fields = map {$_ => 1} ref($config->{fields}) eq 'ARRAY' ? @{$config->{fields}} : split(',', $config->{fields});

  # do some trickery to make sure we're not running unnecessary code
  # this first one only switches on the HGVS options for the requested fields  
  $config->{$_} = 1 for grep {$_ =~ /^hgvs/} keys %set_fields;

  # and this one switches on check_existing if the user wants variant IDs
  my %opt_map = ('id' => 'check_existing');
  $config->{$opt_map{$_}} = 1 for grep {$set_fields{$_}} keys %opt_map;

  # set up/down distance to 0, we only want overlaps
  $config->{distance} = 0;
  
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

  while(my $line = $self->next_output_line(1)) {
    delete($line->{id});
    my $line_id = $line->{input};
    
    merge_arrays($order, [$line_id]);
    find_in_ref($line, \%want_keys, $results->{$line_id} ||= {input => $line_id});

    if(@{$self->warnings}) {
      $results->{$line_id}->{warnings} = [map {$_->{msg}} @{$self->warnings}];
      $self->internalise_warnings();
    }
  }

  return [map {$results->{$_}} @$order];
}

1;
