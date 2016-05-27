# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::VEP::Utils qw(format_coords convert_arrayref numberify merge_hashes get_time);

## format_coords
################

is(format_coords(1, 1), '1',   'format_coords - same');
is(format_coords(1, 2), '1-2', 'format_coords - diff 1');
is(format_coords(2, 1), '1-2', 'format_coords - diff 2');

is(format_coords(1), '1-?', 'format_coords - missing 1');
is(format_coords(undef, 1), '?-1', 'format_coords - missing 2');
is(format_coords(undef, undef), '-', 'format_coords - missing 3');


## convert_arrayref
###################

is(convert_arrayref('foo'), 'foo', 'convert_arrayref - scalar');
is(convert_arrayref(['foo', 'bar']), 'foo,bar', 'convert_arrayref - arrayref');
is(convert_arrayref(['foo', 'bar'], '&'), 'foo&bar', 'convert_arrayref - arrayref with separator');


## numberify
############

is_deeply(
  numberify(['0', '1']),
  [0, 1],
  'numberify - arrayref'
);

is_deeply(
  numberify({foo => '0', bar => '0.1'}),
  {foo => 0, bar => 0.1},
  'numberify - hashref'
);

is_deeply(
  numberify({foo => '0', bar => ['1', '2']}),
  {foo => 0, bar => [1,2]},
  'numberify - mixed'
);

is_deeply(
  numberify(['0.0001']),
  [0.0001],
  'numberify - float'
);

is_deeply(
  numberify(['1e3']),
  [1000],
  'numberify - exponential'
);

is_deeply(
  numberify({id => '123', seq_region_name => '123'}, {id => 1, seq_region_name => 1}),
  {id => '123', seq_region_name => '123'},
  'numberify - exempt keys intact'
);



## merge_hashes
###############

is_deeply(
  merge_hashes({a => 1}, {b => 2}),
  {a => 1, b => 2},
  'merge_hashes - simple'
);

is_deeply(
  merge_hashes({a => 1, b => {c => [1, 2], d => 3}}, {b => {e => 5}, f => 6}),
  {a => 1, b => {c => [1, 2], d => 3, e => 5}, f => 6},
  'merge_hashes - complex'
);

is_deeply(
  merge_hashes({a => 1}, {a => 2}),
  {a => 2},
  'merge_hashes - overwrite'
);

is_deeply(
  merge_hashes({a => [1, 2]}, {a => [3, 4]}),
  {a => [1, 2, 3, 4]},
  'merge_hashes - array'
);


## get_time
###########

ok(get_time =~ /\d{4}(\-\d\d){2} \d\d(\:\d\d){2}/, 'get_time');

# done
done_testing();
