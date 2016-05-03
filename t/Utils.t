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

use Bio::EnsEMBL::VEP::Utils qw(format_coords convert_arrayref get_time);

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


## get_time
###########

ok(get_time =~ /\d{4}(\-\d\d){2} \d\d(\:\d\d){2}/, 'get_time');

# done
done_testing();
