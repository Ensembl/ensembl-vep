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
use FindBin qw($Bin);
use Bio::EnsEMBL::VEP::Utils qw(get_version_data);

# set some vars we'll use
my $vep_module   = 'ensembl-vep';
my $version_dir  = $Bin.'/../.version/';
my $version_file = $version_dir.$vep_module;

SKIP: {
  ## REMEMBER TO UPDATE THIS SKIP NUMBER IF YOU ADD MORE TESTS!!!!
  skip 'Not a git repository', 3 unless -d "$Bin/../.git";

  # find out what branch we are on
  my $git_branch = `cd $Bin; git rev-parse --abbrev-ref HEAD`;
  chomp($git_branch);
  ok($git_branch, 'get current git branch');

  # no need to continue unless we're on a release/NN branch
  my $git_release_number;

  if($git_branch !~ m/^release\/(\d+)/) {
    done_testing();
    exit(0);
  }
  else {
    $git_release_number = $1;
  }

  my $software_version_data = get_version_data($version_dir);

  is(ref($software_version_data->{$vep_module}), 'HASH', $vep_module.' version data');
  is($software_version_data->{$vep_module}->{release}, $git_release_number, $vep_module.' release matches git release');
}


done_testing();
