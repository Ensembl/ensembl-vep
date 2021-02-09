# Copyright [2016-2021] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Exception;
use FindBin qw($Bin);

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::FilterSet');

my $fs = Bio::EnsEMBL::VEP::FilterSet->new();
is(ref($fs), 'Bio::EnsEMBL::VEP::FilterSet', 'new ref');



## METHOD TESTS
###############

# defined_and_non_empty
is($fs->defined_and_non_empty(1), 1, 'defined_and_non_empty - 1');
is($fs->defined_and_non_empty(0), 1, 'defined_and_non_empty - 0');
is($fs->defined_and_non_empty('foo'), 1, 'defined_and_non_empty - string');
is($fs->defined_and_non_empty(""), 0, 'defined_and_non_empty - empty string');
is($fs->defined_and_non_empty(undef), 0, 'defined_and_non_empty - undef');
is($fs->defined_and_non_empty(), 0, 'defined_and_non_empty - no value');


# create_filter_node
is_deeply($fs->create_filter_node(), {logic => 'and', components => []}, 'create_filter_node - no params');
is_deeply($fs->create_filter_node({logic => 'or'}), {logic => 'or', components => []}, 'create_filter_node - logic');
is_deeply($fs->create_filter_node({components => ['foo']}), {logic => 'and', components => ['foo']}, 'create_filter_node - components');


# finish_filter_node
is_deeply(
  $fs->finish_filter_node($fs->create_filter_node),
  { logic => 'and', components => [], operator => 'ex' },
  'finish_filter_node - basic'
);

is_deeply(
  $fs->finish_filter_node($fs->create_filter_node({operator => 'gt'})),
  { logic => 'and', components => [], operator => 'gt' },
  'finish_filter_node - dont overwrite op'
);

is_deeply(
  $fs->finish_filter_node($fs->create_filter_node({operator => 'is'})),
  { logic => 'and', components => [], operator => 'ex' },
  'finish_filter_node - no value changes op 1'
);

is_deeply(
  $fs->finish_filter_node($fs->create_filter_node({operator => 'ne'})),
  { logic => 'and', components => [], operator => 'nex' },
  'finish_filter_node - no value changes op 2'
);

is_deeply(
  $fs->finish_filter_node($fs->create_filter_node({parent => 'foo'})),
  { logic => 'and', components => [], operator => 'ex' },
  'finish_filter_node - parent gets deleted'
);

throws_ok
  { $fs->finish_filter_node($fs->create_filter_node({operator => 'in'})) }
  qr/No list\/file given/,
  'finish_filter_node - in operator throws on no value';

throws_ok
  { $fs->finish_filter_node($fs->create_filter_node({operator => 'in', value => 'foo'})) }
  qr/Could not find\/parse list /,
  'finish_filter_node - in operator throws on invalid value';

is_deeply(
  $fs->finish_filter_node($fs->create_filter_node({operator => 'in', value => 'foo,bar'})),
  { logic => 'and', components => [], operator => 'in', value => {foo => 1, bar => 1} },
  'finish_filter_node - in operator with list'
);

is_deeply(
  $fs->finish_filter_node($fs->create_filter_node({operator => 'in', value => $test_cfg->{filter_list}})),
  { logic => 'and', components => [], operator => 'in', value => {'MRPL39' => 1, 'JAM2' => 1} },
  'finish_filter_node - in operator with file'
);


# parse_filters
is_deeply($fs->parse_filters(), {is_root => 1, components => [], logic => 'and'}, 'parse_filters - empty');

is_deeply(
  $fs->parse_filters(['foo']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' }
    ],
    logic => 'and'
  },
  'parse_filters - basic'
);

is_deeply(
  $fs->parse_filters(['not foo']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and', not => 1 }
    ],
    logic => 'and'
  },
  'parse_filters - basic not'
);

is_deeply(
  $fs->parse_filters(['foo', 'bar']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' },
      { components => [], operator => 'ex', field => 'bar', logic => 'and' }
    ],
    logic => 'and'
  },
  'parse_filters - basic, two filters'
);

is_deeply(
  $fs->parse_filters(['foo eq bar'])->{components}->[0],
  { components => [], operator => 'eq', field => 'foo', logic => 'and', value => 'bar' },
  'parse_filters - basic with op and value'
);

is_deeply(
  $fs->parse_filters(['foo = bar'])->{components}->[0],
  { components => [], operator => 'eq', field => 'foo', logic => 'and', value => 'bar' },
  'parse_filters - op synonym 1'
);

is_deeply(
  $fs->parse_filters(['foo is bar'])->{components}->[0],
  { components => [], operator => 'eq', field => 'foo', logic => 'and', value => 'bar' },
  'parse_filters - op synonym 2'
);

is_deeply(
  $fs->parse_filters(['foo and bar']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' },
      { components => [], operator => 'ex', field => 'bar', logic => 'and' }
    ],
    logic => 'and'
  },
  'parse_filters - two with logic (and)'
);

is_deeply(
  $fs->parse_filters(['foo or bar']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' },
      { components => [], operator => 'ex', field => 'bar', logic => 'or' }
    ],
    logic => 'and'
  },
  'parse_filters - two with logic (or)'
);

is_deeply(
  $fs->parse_filters(['foo eq goo and bar']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'eq', field => 'foo', logic => 'and', value => 'goo' },
      { components => [], operator => 'ex', field => 'bar', logic => 'and' }
    ],
    logic => 'and'
  },
  'parse_filters - two, one short one full 1'
);

is_deeply(
  $fs->parse_filters(['foo and bar eq goo']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' },
      { components => [], operator => 'eq', field => 'bar', logic => 'and', value => 'goo' }
    ],
    logic => 'and'
  },
  'parse_filters - two, one short one full 2'
);

is_deeply(
  $fs->parse_filters(['foo and (bar)']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' },
      { operator => 'ex', logic => 'and', components => [
        { components => [], operator => 'ex', field => 'bar', logic => 'and' }
      ] },
    ],
    logic => 'and'
  },
  'parse_filters - parentheses basic'
);

is_deeply(
  $fs->parse_filters(['foo and not (bar or gwar)']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' },
      { operator => 'ex', logic => 'and', not => 1, components => [
        { components => [], operator => 'ex', field => 'bar', logic => 'and' },      
        { components => [], operator => 'ex', field => 'gwar', logic => 'or' }
      ] },
    ],
    logic => 'and'
  },
  'parse_filters - parentheses basic with not'
);

is_deeply(
  $fs->parse_filters(['foo and (bar or (gwar and char)) or (car and star)']),
  {
    is_root => 1,
    components => [
      { components => [], operator => 'ex', field => 'foo', logic => 'and' },
      { operator => 'ex', logic => 'and', components => [
        { components => [], operator => 'ex', field => 'bar', logic => 'and' },      
        { operator => 'ex', logic => 'or', components => [
          { components => [], operator => 'ex', field => 'gwar', logic => 'and' },
          { components => [], operator => 'ex', field => 'char', logic => 'and' },
        ] },
      ] },
      { operator => 'ex', logic => 'or', components => [
        { components => [], operator => 'ex', field => 'car', logic => 'and' },      
        { components => [], operator => 'ex', field => 'star', logic => 'and' }
      ] },
    ],
    logic => 'and'
  },
  'parse_filters - parentheses complex'
);

throws_ok
  { $fs->parse_filters(['(']) }
  qr/Error parsing filter string.+incomplete parentheses sets/,
  'parse_filters - incomplete parentheses 1';

throws_ok
  { $fs->parse_filters(['foo and (bar']) }
  qr/Error parsing filter string.+incomplete parentheses sets/,
  'parse_filters - incomplete parentheses 2';

throws_ok
  { $fs->parse_filters([')']) }
  qr/Error parsing filter string.+incomplete parentheses sets/,
  'parse_filters - incomplete parentheses 3';

throws_ok
  { $fs->parse_filters(['foo and (bar or (gwar and char) or (car and star)']) }
  qr/Error parsing filter string.+incomplete parentheses sets/,
  'parse_filters - incomplete parentheses 4';

throws_ok
  { $fs->parse_filters(['foo was bar']) }
  qr/No such operator/,
  'parse_filters - invalid operator';



# new
is_deeply(
  Bio::EnsEMBL::VEP::FilterSet->new('foo'),
  bless( {
    'filter_root' => {
      'components' => [
        { components => [], operator => 'ex', field => 'foo', logic => 'and' }
      ],
      'logic' => 'and',
      'is_root' => 1
    }
  }, 'Bio::EnsEMBL::VEP::FilterSet' ),
  'new with filter'
);


# evaluate - ex
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo');
is($fs->evaluate({foo => 1}),     1, 'evaluate - ex - pass');
is($fs->evaluate({foo => 0}),     1, 'evaluate - ex - pass 0');
is($fs->evaluate({bar => 1}),     0, 'evaluate - ex - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - ex - fail undef');

# evaluate - nex
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo nex');
is($fs->evaluate({bar => 1}),     1, 'evaluate - nex - pass');
is($fs->evaluate({foo => 0}),     0, 'evaluate - nex - fail 0');
is($fs->evaluate({foo => undef}), 1, 'evaluate - nex - pass undef');

# evaluate - eq
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo eq 1');
is($fs->evaluate({foo => 1}),     1, 'evaluate - eq - pass');
is($fs->evaluate({foo => 11}),    0, 'evaluate - eq - fail with containing');
is($fs->evaluate({foo => 0}),     0, 'evaluate - eq - fail 0');
is($fs->evaluate({bar => 1}),     0, 'evaluate - eq - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - eq - fail undef');

# evaluate - ne
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo ne 1');
is($fs->evaluate({foo => 0}),     1, 'evaluate - ne - pass');
is($fs->evaluate({foo => 1}),     0, 'evaluate - ne - fail 0');
is($fs->evaluate({bar => 1}),     0, 'evaluate - ne - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - ne - fail undef');

# evaluate - gt
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo gt 5');
is($fs->evaluate({foo => 10}),    1, 'evaluate - gt - pass');
is($fs->evaluate({foo => 0}),     0, 'evaluate - gt - fail 0');
is($fs->evaluate({foo => 5}),     0, 'evaluate - gt - fail match');
is($fs->evaluate({foo => -1}),    0, 'evaluate - gt - fail negative');
is($fs->evaluate({foo => 1e-5}),  0, 'evaluate - gt - fail scientific');
is($fs->evaluate({foo => 1e5}),   1, 'evaluate - gt - pass scientific');
is($fs->evaluate({foo => 0.1}),   0, 'evaluate - gt - fail decimal');
is($fs->evaluate({foo => 5.1}),   1, 'evaluate - gt - pass decimal');
is($fs->evaluate({bar => 1}),     0, 'evaluate - gt - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - gt - fail undef');
is($fs->evaluate({foo => '10&15'}), 1, 'evaluate - gt - pass multiple values');
is($fs->evaluate({foo => '1&4'}), 0, 'evaluate - gt - fail multiple values');

# evaluate - gte
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo gte 5');
is($fs->evaluate({foo => 10}),    1, 'evaluate - gte - pass');
is($fs->evaluate({foo => 0}),     0, 'evaluate - gte - fail 0');
is($fs->evaluate({foo => 5}),     1, 'evaluate - gte - pass match');
is($fs->evaluate({foo => -1}),    0, 'evaluate - gte - fail negative');
is($fs->evaluate({foo => 1e-5}),  0, 'evaluate - gte - fail scientific');
is($fs->evaluate({foo => 1e5}),   1, 'evaluate - gte - pass scientific');
is($fs->evaluate({foo => 0.1}),   0, 'evaluate - gte - fail decimal');
is($fs->evaluate({foo => 5.1}),   1, 'evaluate - gte - pass decimal');
is($fs->evaluate({bar => 1}),     0, 'evaluate - gte - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - gte - fail undef');

# evaluate - lt
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo lt 5');
is($fs->evaluate({foo => 1}),     1, 'evaluate - lt - pass');
is($fs->evaluate({foo => 0}),     1, 'evaluate - lt - pass 0');
is($fs->evaluate({foo => 5}),     0, 'evaluate - lt - fail match');
is($fs->evaluate({foo => -1}),    1, 'evaluate - lt - pass negative');
is($fs->evaluate({foo => 1e5}),   0, 'evaluate - lt - fail scientific');
is($fs->evaluate({foo => 1e-5}),  1, 'evaluate - lt - pass scientific');
is($fs->evaluate({foo => 5.1}),   0, 'evaluate - lt - fail decimal');
is($fs->evaluate({foo => 0.1}),   1, 'evaluate - lt - pass decimal');
is($fs->evaluate({bar => 1}),     0, 'evaluate - lt - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - lt - fail undef');

# evaluate - lte
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo lte 5');
is($fs->evaluate({foo => 1}),     1, 'evaluate - lte - pass');
is($fs->evaluate({foo => 0}),     1, 'evaluate - lte - pass 0');
is($fs->evaluate({foo => 5}),     1, 'evaluate - lte - pass match');
is($fs->evaluate({foo => -1}),    1, 'evaluate - lte - pass negative');
is($fs->evaluate({foo => 1e5}),   0, 'evaluate - lte - fail scientific');
is($fs->evaluate({foo => 1e-5}),  1, 'evaluate - lte - pass scientific');
is($fs->evaluate({foo => 5.1}),   0, 'evaluate - lte - fail decimal');
is($fs->evaluate({foo => 0.1}),   1, 'evaluate - lte - pass decimal');
is($fs->evaluate({bar => 1}),     0, 'evaluate - lte - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - lte - fail undef');

# evaluate - re
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo re bar');
is($fs->evaluate({foo => 'bar'}),  1, 'evaluate - re - pass match');
is($fs->evaluate({foo => 'bare'}), 1, 'evaluate - re - pass part match 1');
is($fs->evaluate({foo => 'ebar'}), 1, 'evaluate - re - pass part match 2');
is($fs->evaluate({foo => 'bear&bar'}), 1, 'evaluate - re - pass with & character');
is($fs->evaluate({foo => 'bear'}), 0, 'evaluate - re - fail');
is($fs->evaluate({bar => 1}),      0, 'evaluate - re - fail missing');
is($fs->evaluate({foo => undef}),  0, 'evaluate - re - fail undef');

$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo re b.r');
is($fs->evaluate({foo => 'bar'}),  1, 'evaluate - re - . wildcard - pass');
is($fs->evaluate({foo => 'br'}),   0, 'evaluate - re - . wildcard - fail no char');
is($fs->evaluate({foo => 'baar'}), 0, 'evaluate - re - . wildcard - fail >1 char');

$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo re b.*r');
is($fs->evaluate({foo => 'bar'}),  1, 'evaluate - re - .* wildcard - pass');
is($fs->evaluate({foo => 'br'}),   1, 'evaluate - re - .* wildcard - pass no char');
is($fs->evaluate({foo => 'baar'}), 1, 'evaluate - re - .* wildcard - pass >1 char');

$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo re b.+r');
is($fs->evaluate({foo => 'bar'}),  1, 'evaluate - re - .+ wildcard - pass');
is($fs->evaluate({foo => 'br'}),   0, 'evaluate - re - .+ wildcard - fail no char');
is($fs->evaluate({foo => 'baar'}), 1, 'evaluate - re - .+ wildcard - pass >1 char');

$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo re b[abcde]r');
is($fs->evaluate({foo => 'bar'}), 1, 'evaluate - re - charset - pass');
is($fs->evaluate({foo => 'br'}),  0, 'evaluate - re - charset - fail no char');
is($fs->evaluate({foo => 'byr'}), 0, 'evaluate - re - charset - fail char not in set');
is($fs->evaluate({foo => 'b1r'}), 0, 'evaluate - re - charset - fail num');

# evaluate - nre
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo nre bar');
is($fs->evaluate({foo => 'car'}),  1, 'evaluate - nre - pass');
is($fs->evaluate({foo => 'bare'}), 0, 'evaluate - nre - fail part match 1');
is($fs->evaluate({foo => 'ebar'}), 0, 'evaluate - nre - fail part match 2');
is($fs->evaluate({foo => 'bear&bar'}), 0, 'evaluate - re - fail with & character');
is($fs->evaluate({foo => 'bear'}), 1, 'evaluate - nre - pass near match');
is($fs->evaluate({bar => 1}),      0, 'evaluate - nre - fail missing');
is($fs->evaluate({foo => undef}),  0, 'evaluate - nre - fail undef');

# evaluate - in
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo in 1,2');
is($fs->evaluate({foo => 1}),     1, 'evaluate - in list - pass 1');
is($fs->evaluate({foo => 2}),     1, 'evaluate - in list - pass 2');
is($fs->evaluate({foo => 3}),     0, 'evaluate - in list - fail');
is($fs->evaluate({bar => 1}),     0, 'evaluate - in list - fail missing');
is($fs->evaluate({foo => undef}), 0, 'evaluate - in list - fail undef');

$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo in '.$test_cfg->{filter_list});
is($fs->evaluate({foo => 'MRPL39'}), 1, 'evaluate - in file - pass 1');
is($fs->evaluate({foo => 'JAM2'}),   1, 'evaluate - in file - pass 2');
is($fs->evaluate({foo => 3}),        0, 'evaluate - in file - fail');
is($fs->evaluate({bar => 1}),        0, 'evaluate - in file - fail missing');
is($fs->evaluate({foo => undef}),    0, 'evaluate - in file - fail undef');

# evaluate - and
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo and bar');
is($fs->evaluate({foo => 1, bar => 1}), 1, 'evaluate - and - pass');
is($fs->evaluate({foo => 1}),           0, 'evaluate - and - fail');

# evaluate - or
$fs = Bio::EnsEMBL::VEP::FilterSet->new('foo or bar');
is($fs->evaluate({foo => 1, bar => 1}), 2, 'evaluate - or - pass 1');
is($fs->evaluate({foo => 1}),           1, 'evaluate - or - pass 2');
is($fs->evaluate({goo => 1}),           0, 'evaluate - or - fail');

# evaluate - not
$fs = Bio::EnsEMBL::VEP::FilterSet->new('not foo');
is($fs->evaluate({bar => 1}),     1, 'evaluate - not - pass');
is($fs->evaluate({foo => undef}), 1, 'evaluate - not - pass undef');
is($fs->evaluate({foo => 1}),     0, 'evaluate - not - fail');

# evaluate - synonyms
is(Bio::EnsEMBL::VEP::FilterSet->new('foo')->evaluate({FOO => 1}),  1, 'evaluate - synonym - uc');
is(Bio::EnsEMBL::VEP::FilterSet->new('FOO')->evaluate({foo => 1}),  1, 'evaluate - synonym - lc');
is(Bio::EnsEMBL::VEP::FilterSet->new('foo')->evaluate({_foo => 1}), 1, 'evaluate - synonym - underscore');

# evaluate - value is field
is(Bio::EnsEMBL::VEP::FilterSet->new('FOO > #BAR')->evaluate({FOO => 2, BAR => 1}),  1, 'evaluate - value is field - pass');
is(Bio::EnsEMBL::VEP::FilterSet->new('FOO > #bar')->evaluate({FOO => 2, BAR => 1}),  1, 'evaluate - value is field - pass synonym');
is(Bio::EnsEMBL::VEP::FilterSet->new('FOO < #BAR')->evaluate({FOO => 2, BAR => 1}),  0, 'evaluate - value is field - fail');
is(Bio::EnsEMBL::VEP::FilterSet->new('FOO = #BAR')->evaluate({FOO => 2}),  0, 'evaluate - value is field - fail missing');


{
  my $fs = Bio::EnsEMBL::VEP::FilterSet->new('foo');

  $fs->evaluate({bar => 1});
  is_deeply($fs->synonyms, {}, 'empty synonyms without limit_synonym_search');

  $fs->limit_synonym_search(1);
  $fs->evaluate({bar => 1});
  is_deeply($fs->synonyms, {foo => undef}, 'undef synonym with limit_synonym_search');
}

## ontology
###########


SKIP: {
  my $db_cfg = $test_cfg->db_cfg;

  eval q{
    use Bio::EnsEMBL::Test::TestUtils;
    use Bio::EnsEMBL::Test::MultiTestDB;
    1;
  };

  my $can_use_db = $db_cfg && scalar keys %$db_cfg && !$@;

  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'No local database configured', 8 unless $can_use_db;

  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('multi') if $can_use_db;

  my $fs = Bio::EnsEMBL::VEP::FilterSet->new('foo is_child gene_variant');

  throws_ok { $fs->evaluate({foo => 1}) } qr/No ontology adaptor set/, 'is_child - no adaptor';

  ok($fs->ontology_adaptor($multi->get_DBAdaptor('ontology')->get_OntologyTermAdaptor), 'is_child - set adaptor');

  throws_ok { $fs->evaluate({foo => 1}) } qr/No ontology_name set/, 'is_child - no ontology_name';

  ok($fs->ontology_name('SO'), 'is_child - set ontology_name');

  is($fs->evaluate({foo => 1}),                  0, 'is_child - fail');
  is($fs->evaluate({foo => 'gene_variant'}),     1, 'is_child - pass exact');
  is($fs->evaluate({foo => 'missense_variant'}), 1, 'is_child - pass child 1');
  is($fs->evaluate({foo => 'intron_variant'}),   1, 'is_child - pass child 2');
}


done_testing();
