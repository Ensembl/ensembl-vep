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
  test_ini_file  => $Bin.'/testdata/vep.ini',
  registry_file  => $Bin.'/testdata/vep.registry',
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

sub registry_file {
  my $self = shift;
  my $dbname = shift;

  if(!exists($self->{user_registry_file})) {

    die("ERROR: No dbname given\n") unless $dbname;

    my %db_cfg = %{$self->db_cfg()};

    die("ERROR: No db config found\n") unless scalar keys %db_cfg;

    $db_cfg{dbname} = $dbname;

    my $base_file = $self->{registry_file};

    # we need to write a new version of this file
    my $reg_file = $base_file.$$;

    open IN, $base_file or die "ERROR: Could not read from $base_file\n";
    open OUT, ">$reg_file" or die "ERROR: Could not write to $reg_file\n";

    while(<IN>) {
      my $line = $_;
      $line =~ s/\_\_$_\_\_/$db_cfg{$_}/ge for keys %db_cfg;
      print OUT $line;
    }

    close IN;
    close OUT;

    $self->{user_registry_file} = $reg_file;
  }

  return $self->{user_registry_file};
}

sub DESTROY {
  my $self = shift;
  unlink($self->{user_registry_file}) if $self->{user_registry_file};
}