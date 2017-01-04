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

# EnsEMBL module for Bio::EnsEMBL::VEP::FilterSet
#
#

=head1 NAME

Bio::EnsEMBL::VEP::FilterSet - Base class used for filters

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::FilterSet;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# package vars
my %FILTER_SYNONYMS = (
  '>' => 'gt',
  '>=' => 'gte',
  '<' => 'lt',
  '<=' => 'lte',
  
  'is' => 'eq',
  '=' => 'eq',
  
  '!=' => 'ne',
  
  'exists' => 'ex',
  'defined' => 'ex',
  
  'match' => 're',
  'matches' => 're',
  'regex' => 're',
);

sub new {  
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # initialise self
  my $self = bless {}, $class;

  $self->{filter_root} = $self->parse_filters(\@_);

  return $self;
}

sub parse_filters {
  my $self = shift;
  my $filters = shift;

  my $root = $self->create_filter_node({
    is_root => 1
  });
  
  foreach my $filter_list(@$filters) {

    my $current = $self->create_filter_node({
      parent => $root,
    });

    while($filter_list =~ m/([^\(^\)^\s]*?)(\s|\(|\)|$)/g) {
      my ($word, $sep) = ($1, $2);

      # no word or separator - should be end of the string
      unless($self->defined_and_non_empty($word) || $self->defined_and_non_empty($sep)) {
        my $parent = $current->{parent};

        # parent should now be root if all opened parentheses have been closed
        throw("ERROR: Error parsing filter string - incomplete parentheses sets?\n") unless $parent->{is_root};

        $self->finish_filter_node($current);

        push @{$parent->{components}}, $current;
      }

      if($self->defined_and_non_empty($word)) {

        if($word eq 'and' || $word eq 'or') {
          my $parent = $current->{parent};

          $self->finish_filter_node($current);
          push @{$parent->{components}}, $current;

          $current = $self->create_filter_node({
            logic => $word,
            parent => $parent,
          });
        }

        # invert with not?
        elsif($word eq 'not') {
          $current->{not} = 1;
        }

        ## current gets field, operator, value in that order

        # process field
        elsif(!$current->{field}) {
          $current->{field} = $word;
        }

        # process operator
        elsif(!$current->{operator}) {

          # check operator is valid
          my $sub_name = 'filter_'.$word;
          
          if(!defined(&$sub_name)) {
            if(defined($FILTER_SYNONYMS{$word})) {
              $word = $FILTER_SYNONYMS{$word};
              $sub_name = 'filter_'.$word;
            }
            else {
              throw("ERROR: No such operator \"$word\"\n") unless defined(&$sub_name);
            }
          }

          $current->{operator} = $word;
        }

        # process value
        elsif(!defined($current->{value})) {
          $current->{value} = $word;
        }
      }

      if($sep) {
        if($sep eq '(') {
          my $parent = $current;

          $current = $self->create_filter_node({
            parent => $parent,
          });
        }
        elsif($sep eq ')') {
          my $parent = $current->{parent};

          # finish child
          $self->finish_filter_node($current);
          push @{$parent->{components}}, $current;
          $current = $parent;
        }
      }
    }
  }
  
  return $root;
}

sub defined_and_non_empty {
  return defined($_[1]) && $_[1] ne '' ? 1 : 0;
}

sub create_filter_node {
  my $self = shift;
  my $filter = shift;

  $filter->{logic} = 'and' unless exists($filter->{logic});
  $filter->{components} = [] unless exists($filter->{components});

  return $filter;
}

sub finish_filter_node {
  my $self = shift;
  my $filter = shift;

  $filter->{operator} ||= 'ex';

  # if(defined($config->{ontology}) && $filter->{field} eq 'Consequence' && $filter->{operator} eq 'eq') {
  #   $filter->{operator} = 'is_child';
  # }

  if(!defined($filter->{value})) {
    $filter->{operator} = 'nex' if $filter->{operator} eq 'ne';
    $filter->{operator} = 'ex' if $filter->{operator} eq 'is';
  }

  if($filter->{operator} eq 'in') {
    my %compare;
    my $list = $filter->{value};

    throw("ERROR: No list/file given\n") unless $list;
  
    if($list =~ /.+\,.+/) {
      %compare = map {$_ => 1} split(',', $list);
    }
    
    # file?
    elsif(-e $list) {
      open IN, $list or throw("ERROR: Could not read from file $list\n");
      while(<IN>) {
        chomp;

        # perl 5.8.8 doesn't recognise \v meaning any vertical whitespace...
        if($] < 5.01) {
          s/\r|(\x0D\x0A).+//;
        }
        else {
          s/\r|(?>\v|\x0D\x0A)//g;
        }

        $compare{$_} = 1;
      }
      close IN;
    }
    
    else {
      throw("ERROR: Could not find/parse list $list\n");
    }

    $filter->{value} = \%compare;
  }

  delete $filter->{parent};

  return $filter;
}

sub evaluate {
  my $self = shift;
  my $data = shift;
  my $filter = @_ ? shift : $self->{filter_root};

  my $return = 1;

  if(scalar @{$filter->{components}}) {

    foreach my $sub(@{$filter->{components}}) {

      my $value = $self->evaluate($data, $sub);

      if($sub->{logic} eq 'and') {
        $return *= $value;
      }
      elsif($sub->{logic} eq 'or') {
        $return += $value;
      } 
    }
  }

  else {
    my $predicate_name = 'filter_'.$filter->{operator};
    my $predicate = \&$predicate_name;

    # process input
    my $field = $filter->{field};
    my $value = $filter->{value};
    my $input = $self->get_input($field, $value, $data);
    
    # run filter
    if(defined($value) && !defined($input)) {
      $return = 0;
    }
    else {
      $return = &$predicate($input, $value, $self);
    }
  }

  if($filter->{not}) {
    $return = $return == 0 ? 1 : 0;
  }

  return $return;
}

sub get_input {
  my ($self, $field, $value, $data) = @_;

  my $input;

  if(exists($data->{$field})) {
    $input = $data->{$field};
  }
  else {

    # try synonyms
    my $synonyms = $self->synonyms;

    if(my $synonym = $synonyms->{$field}) {
      $input = $data->{$synonym};
    }

    else {
      # lc
      if(exists($data->{lc($field)})) {
        $synonyms->{$field} = lc($field);
        $input = $data->{$synonyms->{$field}};
      }
      # uc
      elsif(exists($data->{uc($field)})) {
        $synonyms->{$field} = uc($field);
        $input = $data->{$synonyms->{$field}};
      }
      # prefix with "_"
      elsif(exists($data->{'_'.$field})) {
        $synonyms->{$field} = '_'.$field;
        $input = $data->{$synonyms->{$field}};
      }
      # substring
      else {
        my @matches = grep {$_ =~ /^$field/i} keys %$data;

        if(scalar @matches == 1) {
          $synonyms->{$field} = $matches[0];
        }
        elsif(scalar @matches > 1) {
          throw(
            sprintf(
              'ERROR: Multiple data fields found matching query field "%s": %s',
              $field,
              join(", ", @matches)
            )
          );
        }
      }
    }
  }

  if(defined($input) && $input =~ /([\w\.\-]+)?\:?\(?([\-\d\.e]*)\)?/ && $field ne 'CELL_TYPE') {

    my ($text, $num) = ($1, $2);
    
    if($value && $value =~ /^[\-\d\.e]+$/) {
      $input = $text =~ /^\-?\d+\.?\d*(e\-?\d+)?$/ ? $text : $num;
    }
    else {
      $input = $text;
    }
  }

  return $input;
}

sub synonyms {
  my $self = shift;
  $self->{synonyms} = shift if @_;
  return $self->{synonyms} ||= {};
}

sub ontology_adaptor {
  my $self = shift;
  $self->{_ontology_adaptor} = shift if @_;
  return $self->{_ontology_adaptor};
}

sub ontology_name {
  my $self = shift;
  $self->{_ontology_name} = shift if @_;
  return $self->{_ontology_name};
}



## PREDICATE METHODS
####################

# basic filters
sub filter_eq  { return $_[0] eq $_[1] }
sub filter_ne  { return $_[0] ne $_[1] }
sub filter_gt  { return $_[0] >  $_[1] }
sub filter_lt  { return $_[0] <  $_[1] }
sub filter_gte { return $_[0] >= $_[1] }
sub filter_lte { return $_[0] <= $_[1] }
sub filter_ex  { return defined($_[0]) }
sub filter_nex { return !defined($_[0]) }

# string
sub filter_re  { return $_[0] =~ /$_[1]/i }
sub filter_nre { return $_[0] !~ /$_[1]/i }

# in list
sub filter_in  { return exists($_[1]->{$_[0]}) }

# use ontology to check if input is a child of value
# requires you to have set ontology_adaptor and ontology_name
sub filter_is_child {
  my ($child, $parent, $obj) = @_;
  
  # exact match, don't need to use ontology
  return 1 if filter_re($child, $parent);

  my $cache = $obj->{_cache}->{descendants} ||= {};
  
  # get parent term and descendants and cache it on $config
  unless($cache->{$parent}) {

    my $oa = $obj->ontology_adaptor;
    my $ontology_name = $obj->ontology_name;

    throw("ERROR: No ontology adaptor set") unless $oa;
    throw("ERROR: No ontology_name set") unless $ontology_name;
    
    # connect to DBs here
    my $terms = $oa->fetch_all_by_name($parent, $ontology_name);
    throw("ERROR: No matching SO terms found for $parent\n") unless $terms && scalar @$terms;
    throw("ERROR: Found more than one SO term matching $parent: ".join(", ", map {$_->name} @$terms)."\n") if scalar @$terms > 1;
    my $parent_term = $terms->[0];
    $cache->{$parent} = [map {$_->name} @{$parent_term->descendants}];
  }  
  
  return grep {lc($_) eq lc($child)} @{$cache->{$parent}};
}

1;