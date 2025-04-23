# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

use File::Basename;
use File::Spec;
use Test::More;

use Bio::EnsEMBL::Test::MultiTestDB;

diag 'Starting database and files cleaning up...';

my $curr_file = __FILE__;
my $db_conf = Bio::EnsEMBL::Test::MultiTestDB->get_db_conf(dirname(__FILE__));

foreach my $species ( keys %{ $db_conf->{'databases'} } ) {
  my $multi = Bio::EnsEMBL::Test::MultiTestDB->new($species);
}

note "Deleting $curr_file";
my $result = unlink $curr_file;
ok($result, 'Unlink of '.$curr_file.' worked');

done_testing();
