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

## BASIC TESTS
##############

# use test
use_ok('Bio::EnsEMBL::VEP::BaseVEP');

my $bv = Bio::EnsEMBL::VEP::BaseVEP->new();
ok($bv, 'new is defined');

is(ref($bv), 'Bio::EnsEMBL::VEP::BaseVEP', 'check class');

throws_ok { Bio::EnsEMBL::VEP::BaseVEP->new('not a hashref') } qr/Can\'t use .+ as a HASH ref/, 'new with invalid arg';

throws_ok { Bio::EnsEMBL::VEP::BaseVEP->new({config => {}}) } qr/was expected to be .*Bio::EnsEMBL::VEP::Config/, 'new with invalid config object';



## METHOD TESTS
###############

# need to get a config object for further tests
use_ok('Bio::EnsEMBL::VEP::Config');

my $cfg = Bio::EnsEMBL::VEP::Config->new();
ok($cfg, 'get new config object');

$bv = Bio::EnsEMBL::VEP::BaseVEP->new({config => $cfg});
ok($bv, 'new with config object');

is(ref($bv->config), 'Bio::EnsEMBL::VEP::Config', 'config method');

is($bv->param('species'), 'homo_sapiens', 'param get');
is($bv->param('species', 'mouse'), 'mouse', 'param set');
throws_ok { $bv->param() } qr/No param/, 'param without key';



## status_msg tests require we mess with STDOUT
###############################################

# status_msg prints to STDOUT
no warnings 'once';
open(SAVE, ">&STDOUT") or die "Can't save STDOUT\n"; 

close STDOUT;
my $tmp;
open STDOUT, '>', \$tmp;

$bv->status_msg('Hello');
ok($tmp =~ /Hello/, 'status_msg');

open(STDOUT, ">&SAVE") or die "Can't restore STDOUT\n";


## DONE
#######
done_testing();

