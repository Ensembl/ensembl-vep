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

use Bio::EnsEMBL::VEP::Utils qw(
  format_coords
  convert_arrayref
  numberify
  merge_hashes
  merge_arrays
  find_in_ref
  get_time
  get_compressed_filehandle
  get_version_data
  get_version_string
);

use FindBin qw($Bin);
use Scalar::Util qw(looks_like_number);
use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

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
  merge_hashes({b => {e => 5}, f => 6}, {a => 1, b => {c => [1, 2], d => 3}}),
  {a => 1, b => {c => [1, 2], d => 3, e => 5}, f => 6},
  'merge_hashes - complex inverted'
);

is_deeply(
  merge_hashes({a => 1}, {a => 2}),
  {a => 2},
  'merge_hashes - overwrite'
);

is_deeply(
  merge_hashes({a => 1}, {a => 2}, 1),
  {a => 3},
  'merge_hashes - add'
);

is_deeply(
  merge_hashes({a => [1, 2]}, {a => [3, 4]}),
  {a => [1, 2, 3, 4]},
  'merge_hashes - array'
);



## merge_arrays
###############

is_deeply(
  merge_arrays([1], [2]),
  [1, 2],
  'merge_arrays - simple'
);

is_deeply(
  merge_arrays([1], [1, 2]),
  [1, 2],
  'merge_arrays - duplicates'
);

is_deeply(
  merge_arrays([1, 3], [2, 1, 3, 4]),
  [1, 3, 2, 4],
  'merge_arrays - order'
);


## find_in_ref
##############

my $want_keys = { foo => 1, bar => 1};

is_deeply(
  find_in_ref(
    { foo => 'hello' },
    $want_keys
  ),
  { foo => ['hello'] },
  'find_in_ref - simple'
);

is_deeply(
  find_in_ref(
    { foo => 'hello', bar => 'bye' },
    $want_keys
  ),
  { foo => ['hello'], bar => ['bye'] },
  'find_in_ref - two keys'
);

is_deeply(
  find_in_ref(
    [ {foo => 'hello'}, {foo => 'bye'} ],
    $want_keys
  ),
  { foo => ['hello', 'bye'] },
  'find_in_ref - multiple values'
);

is_deeply(
  find_in_ref(
    { foo => ['hello', 'bye'] },
    $want_keys
  ),
  { foo => ['hello', 'bye'] },
  'find_in_ref - multiple values in arrayref'
);

is_deeply(
  find_in_ref(
    { foo => [ {foo => 'hello'}, {bar => 'bye'} ] },
    $want_keys
  ),
  { foo => ['hello'], bar => ['bye'] },
  'find_in_ref - dont get non-scalars'
);


## get_compressed_filehandle
############################

ok(get_compressed_filehandle($test_cfg->{test_gzvcf}), 'get_compressed_filehandle');
throws_ok {get_compressed_filehandle()} qr/No file/, 'get_compressed_filehandle - no file';
throws_ok {get_compressed_filehandle('foobargoobar')} qr/File .+ does not exist/, 'get_compressed_filehandle - file does not exist';
throws_ok {get_compressed_filehandle($test_cfg->{test_vcf})} qr/File .+ binary/, 'get_compressed_filehandle - file not compressed';

my @bak = (
  $Bio::EnsEMBL::VEP::Utils::CAN_USE_PERLIO_GZIP,
  $Bio::EnsEMBL::VEP::Utils::CAN_USE_IO_UNCOMPRESS,
  $Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP
);

SKIP: {
  skip 'PerlIO::Gzip or gzip not installed', 3 unless $Bio::EnsEMBL::VEP::Utils::CAN_USE_PERLIO_GZIP && $Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP;

  is(
    ref(get_compressed_filehandle($test_cfg->{test_gzvcf})),
    'GLOB',
    'get_compressed_filehandle - PerlIO::Gzip'
  );
  is(
    ref(get_compressed_filehandle($test_cfg->{test_gzvcf}, 1)),
    'GLOB',
    'get_compressed_filehandle - PerlIO::Gzip with multistream uses gzip'
  );

  $Bio::EnsEMBL::VEP::Utils::CAN_USE_PERLIO_GZIP = 0;
  is(
    ref(get_compressed_filehandle($test_cfg->{test_gzvcf})),
    'GLOB',
    'get_compressed_filehandle - gzip'
  );
}

SKIP: {
  skip 'gzip or IO::Uncompress::Gunzip not installed', 1 unless $Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP && $Bio::EnsEMBL::VEP::Utils::CAN_USE_IO_UNCOMPRESS;

  $Bio::EnsEMBL::VEP::Utils::CAN_USE_PERLIO_GZIP = 0;
  $Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP = 0;
  is(
    ref(get_compressed_filehandle($test_cfg->{test_gzvcf})),
    'IO::Uncompress::Gunzip',
    'get_compressed_filehandle - IO::Uncompress::Gunzip'
  );
}

$Bio::EnsEMBL::VEP::Utils::CAN_USE_PERLIO_GZIP = 0;
$Bio::EnsEMBL::VEP::Utils::CAN_USE_IO_UNCOMPRESS = 0;
$Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP = 0;

throws_ok {get_compressed_filehandle($test_cfg->{test_gzvcf})} qr/Cannot read from compressed or binary file/, 'get_compressed_filehandle - no available library';

(
  $Bio::EnsEMBL::VEP::Utils::CAN_USE_PERLIO_GZIP,
  $Bio::EnsEMBL::VEP::Utils::CAN_USE_IO_UNCOMPRESS,
  $Bio::EnsEMBL::VEP::Utils::CAN_USE_GZIP
) = @bak;




## get_time
###########

ok(get_time =~ /\d{4}(\-\d\d){2} \d\d(\:\d\d){2}/, 'get_time');



## version data
###############

my $vd = get_version_data($Bin.'/testdata/version');
my $vep_vd = delete($vd->{'ensembl-vep'});

is(ref($vep_vd), 'HASH', 'get_version_data - ensembl-vep hashref');
ok(looks_like_number($vep_vd->{release}), 'get_version_data - ensembl-vep release');
ok(looks_like_number($vep_vd->{sub}), 'get_version_data - ensembl-vep sub');

is_deeply(
  $vd,
  {
    'ensembl-variation' => {
      'sub' => '0000000',
      'release' => '86'
    },
    'ensembl' => {
      'sub' => '0000001',
      'release' => '86'
    }
  },
  'get_version_data - the rest'
);

ok(
  get_version_string($Bin.'/testdata/version') =~
    /ensembl\s+\: 86\.0000001\s+ensembl\-variation\s+\: 86\.0000000\s+ensembl\-vep\s+\: \d+\.\d+/,
  'get_version_string'
);

# done
done_testing();
