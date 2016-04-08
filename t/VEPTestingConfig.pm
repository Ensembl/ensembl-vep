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

package VEPTestingConfig;

use FindBin qw($Bin);

our %DEFAULTS = (
  cache_root_dir => $Bin.'/testdata/cache/',
  cache_species  => 'homo_sapiens',
  cache_version  => 84,
  cache_assembly => 'GRCh38',
  cache_dir      => $Bin.'/testdata/cache/homo_sapiens/84_GRCh38',
  test_ini_file  => $Bin.'/testdata/test_vep.ini',
);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my %config = %DEFAULTS;

  # initialise self
  my $self = bless \%config, $class;

  return $self;
}

sub base_testing_cfg {
  my $self = shift;

  return {
    dir           => $self->{cache_root_dir},
    species       => $self->{cache_species},
    cache_version => $self->{cache_version},
    assembly      => $self->{cache_assembly},
    offline       => 1,
  }
}

sub db_cfg {
  my $self = shift;

  if(!exists($self->{db_cfg})) {

    my $cfg = {};

    if(open IN, $Bin.'/MultiTestDB.conf') {
      my @lines = <IN>;
      $cfg = eval join('', @lines);
      $cfg->{password} = $cfg->{pass};
      close IN;
    }

    $self->{db_cfg} = $cfg;
  }

  return $self->{db_cfg};
}