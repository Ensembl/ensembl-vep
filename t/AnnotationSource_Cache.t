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
use FindBin qw($Bin);

use lib $Bin;
use VEPTestingConfig;
my $test_cfg = VEPTestingConfig->new();

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::AnnotationSource::Cache');

my $root_dir = $test_cfg->{cache_root_dir};

my $c = Bio::EnsEMBL::VEP::AnnotationSource::Cache->new({root_dir => $root_dir});
ok($c, 'new is defined');


## METHODS
##########

is($c->serializer_type, 'storable', 'serializer_type');
is($c->file_suffix, 'gz', 'file_suffix');

is($c->root_dir, $root_dir, 'root_dir get');
is($c->root_dir('/tmp'), '/tmp', 'root_dir set');
$c->root_dir($root_dir);

done_testing();
