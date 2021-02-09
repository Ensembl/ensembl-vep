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

# EnsEMBL module for Bio::EnsEMBL::VEP::FilterSet
#
#

=head1 NAME

Bio::EnsEMBL::VEP::FilterSet - Base class used for filters

=head1 SYNOPSIS

my $fs = Bio::EnsEMBL::VEP::FilterSet->new("field1 = foo and field2 < 10");

my $pass = $fs->evaluate({field1 => 'foo', field2 => 6});
my $fail = $fs->evaluate({field1 => 'bar', field2 => 12});

=head1 DESCRIPTION

The FilterSet class allows logical filtering by determining whether a given
hashref of data "passes" a set of filters described in a fairly simple
string format.

The filter string consists of one or more units linked together by logical
conditions (AND or OR) and optionally compartmentalised by parentheses.

A unit consists of three components: a field, an operator and a value

 - field: key name as found in the hashref to be evaluated. If no exact match
   is found, the code attempts "synonym" lookups by looking for
   case-insensitive and/or partial string matches e.g. "foobar" => "FooBar",
   "cons" => "consequence"

 - operator: the operator used to compare the value given in the filter string
   to the user-supplied value for field in the evaluated hashref

 - value: the value used to compare. If prefixed with "#", corresponds to
   the value of #field in the data e.g. field1 < #field2

If operator and value are excluded, then only the existence of field is
checked.

A unit may be "inverted" by adding "not" before the unit.

Available operators:

 - eq:       strings equal (eq in perl)
 - ne:       strings inequal (ne in perl)
 - gt:       numerical greater than (> in perl)
 - lt:       numerical less than (< in perl)
 - gte:      numerical greater than or equal to (>= in perl)
 - lte:      numerical less than or equal to (>= in perl)
 - ex:       value is defined (defined in perl)
 - nex:      value is not defined (!defined in perl)
 - re:       regular expression match
 - nre:      regular expression non-match
 - in:       value exists in comma-separated list or file
 - is_child: uses ontology_adaptor to lookup via ontology

=head1 METHODS

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


=head2 new

  Arg 1      : string $filter_string
  Example    : $fs = Bio::EnsEMBL::VEP::FilterSet->new($filter_string);
  Description: Constructor for Bio::EnsEMBL::VEP::FilterString.
  Returntype : Bio::EnsEMBL::VEP::FilterString
  Exceptions : none
  Caller     : filter_vep, AnnotationSource
  Status     : Stable

=cut

sub new {  
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  # initialise self
  my $self = bless {}, $class;

  $self->{filter_root} = $self->parse_filters(\@_);

  return $self;
}


=head2 parse_filters

  Arg 1      : string $filter_string
  Example    : $filter_root = $fs->parse_filters($filter_string);
  Description: Takes a filter string and returns a hashref representing
               the filter structure. The structure consists of a series
               of nodes, each having one or more child components,
               where siblings are logically linked filter units within
               the same level of parentheses.

               Examples:

               "foo is bar" =>
               {
                 is_root => 1,
                 logic => 'and',
                 components => [
                   {
                     logic => 'and',
                     field => 'foo',
                     operator => 'eq',
                     value => 'bar'
                   }
                 ]
               }

               "foo is bar and (a < 10 or b > 100)" =>
               {
                 is_root => 1,
                 logic => 'and',
                 components => [
                   {
                     logic => 'and',
                     field => 'foo',
                     operator => 'eq',
                     value => 'bar',
                   },
                   {
                     logic => 'and',
                     operator => 'ex',
                     components => [
                       {
                         logic => 'and',
                         field => 'a',
                         operator => 'lt',
                         value => 10
                       },
                       {
                         logic => 'or',
                         field => 'b',
                         operator => 'gt',
                         value => 100
                       },
                     ],
                   }
                 ]
               }
  Returntype : hashref
  Exceptions : throws on failure to parse filter string
  Caller     : filter_vep, AnnotationSource
  Status     : Stable

=cut

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


=head2 defined_and_non_empty

  Arg 1      : scalar $value
  Examples   : $ok = $fs->defined_and_non_empty(1)
               $ok = $fs->defined_and_non_empty(0)
               $not_ok = $fs->defined_and_non_empty('')
               $not_ok = $fs->defined_and_non_empty()
  Description: Checks if a given value is defined and not an empty string
  Returntype : bool
  Exceptions : none
  Caller     : parse_filters()
  Status     : Stable

=cut

sub defined_and_non_empty {
  return defined($_[1]) && $_[1] ne '' ? 1 : 0;
}


=head2 create_filter_node

  Arg 1      : hashref $filter_node
  Example    : $node = $fs->create_filter_node()
  Description: Initialises a filter node
  Returntype : hashref
  Exceptions : none
  Caller     : parse_filters()
  Status     : Stable

=cut

sub create_filter_node {
  my $self = shift;
  my $filter = shift;

  $filter->{logic} = 'and' unless exists($filter->{logic});
  $filter->{components} = [] unless exists($filter->{components});

  return $filter;
}


=head2 finish_filter_node

  Arg 1      : hashref $filter_node
  Example    : $node = $fs->finish_filter_node()
  Description: Finalises a filter node; assigns default operator (ex),
               checks value or reads file if operator is "in", deletes
               "parent" key to avoid circular references
  Returntype : hashref
  Exceptions : none
  Caller     : parse_filters()
  Status     : Stable

=cut

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


=head2 evaluate

  Arg 1      : hashref $data
  Arg 2      : (optional) hashref $filter_root
  Example    : $pass = $fs->evaluate($data)
  Description: Evaluates a given hashref of data using the filter structure
               defined by this FilterSet. Called iteratively by child nodes.
  Returntype : bool
  Exceptions : none
  Caller     : user
  Status     : Stable

=cut

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
    my $result;

    my $predicate_name = 'filter_'.$filter->{operator};
    my $predicate = \&$predicate_name;

    # process input
    my $field = $filter->{field};
    my $value = $filter->{value};
    my $input = $self->get_input($field, $value, $data);

    # value is field
    if(defined($value) && substr($value, 0, 1) eq '#') {
      $value = $self->get_input(substr($value, 1), undef, $data);
      return 0 unless defined($value);
    }
    
    # run filter
    if(defined($value) && !defined($input)) {
      $return = 0;
    }
    else {

      # Filter multiple values
      if(defined($input) && $input =~ /[0-9](\,|\&)/) {
        my $delimiter = $1;
        my $found = 0;
        my @split_input = split /$delimiter/, $input;
        foreach my $value_i (@split_input) {
          my $predicted = &$predicate($value_i, $value, $self);
          if($predicted) {
            $found = 1;
          }
        }
        $result = $found;
      }
      else {
        $result = &$predicate($input, $value, $self);
      }

      $return = $result;
    }
  }

  if($filter->{not}) {
    $return = $return == 0 ? 1 : 0;
  }

  return $return;
}


=head2 get_input
  
  Arg 1      : string $field
  Arg 2      : scalar $value
  Arg 3      : hashref $data
  Example    : $input = $fs->get_input('foo', 0.1, {foo => 'bar(0.7)'});
  Description: Gets the value for field from given data hashref. The value
               to be compared against is also supplied such that code can
               work out whether to extract numerical or text part of mixed-type
               values (e.g. 0.7 is returned in example above)
  Returntype : scalar
  Exceptions : throws if field name has partial match to multiple fields in data
               during synonym lookup
  Caller     : evaluate()
  Status     : Stable

=cut

sub get_input {
  my ($self, $field, $value, $data) = @_;

  my $input;

  if(exists($data->{$field})) {
    $input = $data->{$field};
  }
  else {

    # try synonyms
    my $synonyms = $self->synonyms;

    if(exists($synonyms->{$field})) {
      $input = defined($synonyms->{$field}) ? $data->{$synonyms->{$field}} : undef;
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
        elsif($self->limit_synonym_search) {
          $synonyms->{$field} = undef;
        }
      }
    }
  }

  if(defined($input) && $input =~ /([\w\.\-]+)?\:?\(?([\-\d\.e]*)\)?/
                     && $field ne 'CELL_TYPE' && $field !~/HGVS/) {

    my ($text, $num) = ($1, $2);
    
    if($num ne '') {
      if($value && $value =~ /^[\-\d\.e]+$/) {
        $input = $text =~ /^\-?\d+\.?\d*(e\-?\d+)?$/ ? $text : $num;
      }
      else {
        $input = $text;
      }
    }
  }

  return $input;
}


=head2 synonyms
  
  Arg 1      : (optional) hashref $synonyms
  Example    : $syns = $fs->synonyms();
  Description: Getter/setter for hashref containing synonyms mapping field
               names as supplied in filter string to those actually found
               in data.
  Returntype : hashref
  Exceptions : none
  Caller     : get_input()
  Status     : Stable

=cut

sub synonyms {
  my $self = shift;
  $self->{synonyms} = shift if @_;
  return $self->{synonyms} ||= {};
}


=head2 limit_synonym_search
  
  Arg 1      : (optional) bool $limit_search
  Example    : $fs->limit_synonym_search(1);
  Description: Getter/setter for setting controlling how synonyms are looked up.
               By enabling this, if a synonym lookup fails the first time, it
               will not be run again. Otherwise, a failed lookup can be repeated
               each time get_input() is called
  Returntype : hashref
  Exceptions : none
  Caller     : get_input()
  Status     : Stable

=cut

sub limit_synonym_search {
  my $self = shift;
  $self->{limit_synonym_search} = shift if @_;
  return $self->{limit_synonym_search} ||= 0;
}


=head2 ontology_adaptor
  
  Arg 1      : (optional) Bio::EnsEMBL::DBSQL::OntologyTermAdaptor $ad
  Example    : $ad = $fs->ontology_adaptor(1);
  Description: Getter/setter for the ontology adaptor to use when looking
               up child terms for the is_child operator
  Returntype : Bio::EnsEMBL::DBSQL::OntologyTermAdaptor
  Exceptions : none
  Caller     : filter_vep, is_child()
  Status     : Stable

=cut

sub ontology_adaptor {
  my $self = shift;
  $self->{_ontology_adaptor} = shift if @_;
  return $self->{_ontology_adaptor};
}


=head2 ontology_name
  
  Arg 1      : (optional) string $name
  Example    : $name = $fs->ontology_name();
  Description: Getter/setter for the ontology name to use when looking
               up child terms for the is_child operator
  Returntype : string
  Exceptions : none
  Caller     : filter_vep, is_child()
  Status     : Stable

=cut

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
