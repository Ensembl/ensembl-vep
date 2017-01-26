=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::InitDump;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use DBI;

my $DEBUG = 0;

sub param_defaults {
  return {
    'refseq' => 0,
    'include_pattern' => '',
    'exclude_pattern' => '',
  };
}

sub fetch_input {
  my $self = shift;
  
  $self->param(
    'species_jobs',
    [
      map {@{$self->get_all_species_jobs_by_server($_)}}
      @{$self->required_param('dump_servers')}
    ]
  );

  return;
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id({}, 1);
  $self->dataflow_output_id($self->param('species_jobs'), 2);
  
  return;
}

sub get_all_species_jobs_by_server {
  my $self = shift;
  my $server = shift;

  my $connection_string = sprintf(
    "DBI:mysql(RaiseError=>1):host=%s;port=%s",
    $server->{host},
    $server->{port}
  );
  
  # connect to DB
  my $dbc = DBI->connect(
    $connection_string, $server->{user}, $server->{pass}
  );
  
  my $version = $self->param('eg_version') || $self->required_param('ensembl_release');

  my $sth = $dbc->prepare(qq{
    SHOW DATABASES LIKE '%\_core\_$version%'
  });
  $sth->execute();
  
  my $db;
  $sth->bind_columns(\$db);
  
  my @dbs;
  push @dbs, $db while $sth->fetch;
  $sth->finish;
  
  # refseq?
  if($self->param('refseq')) {
    $sth = $dbc->prepare(qq{
      SHOW DATABASES LIKE '%\_otherfeatures\_$version%'
    });
    $sth->execute();
    $sth->bind_columns(\$db);
    
    push @dbs, $db while $sth->fetch;
    $sth->finish;
  }

  # remove master and coreexpression
  @dbs = grep {$_ !~ /master|express/} @dbs;

  # filter on pattern if given
  my $pattern = exists($server->{include_pattern}) ? $server->{include_pattern} : $self->param('include_pattern');
  my $exclude = exists($server->{exclude_pattern}) ? $server->{exclude_pattern} : $self->param('exclude_pattern');
  @dbs = grep {$_ =~ /$pattern/i} @dbs if $pattern;
  @dbs = grep {$_ !~ /$exclude/i} @dbs if $exclude;

  my @return;

  foreach my $current_db_name (@dbs) {

    next if $self->is_strain($dbc, $current_db_name) && $current_db_name =~ /mus_musculus/;

    my $group = 'core';
    
    # special case otherfeatures
    if($current_db_name =~ /otherfeatures/) {

      # check it has refseq transcripts
      $sth = $dbc->prepare(qq{
        SELECT COUNT(*)
        FROM $current_db_name\.transcript
        WHERE stable_id LIKE 'NM%'
        OR source = 'refseq'
      });
      $sth->execute;
      
      my $count;
      $sth->bind_columns(\$count);
      $sth->fetch;
      $sth->finish();
      next unless $count;

      $group = 'otherfeatures';
    }

    my $species_ids = $self->get_species_id_hash($dbc, $current_db_name);
    
    # do we have a variation DB?
    my $var_db_name = $self->has_var_db($dbc, $current_db_name);
    
    # do we have a regulation DB?
    my $reg_db_name = $self->has_reg_build($dbc, $current_db_name);

    my $species_count = 0;
    
    foreach my $species_id(keys %$species_ids) {
      my $assembly = $self->get_assembly($dbc, $current_db_name, $species_id);
      next unless $assembly;
      
      # copy server details
      my %species_hash = %$server;
      
      $species_hash{species} = $species_ids->{$species_id};
      $species_hash{species_id} = $species_id;
      $species_hash{assembly} = $assembly;
      $species_hash{dbname} = $current_db_name;
      $species_hash{group} = $group;
      $species_hash{is_multispecies} = scalar keys %$species_ids > 1 ? 1 : 0;
      $species_hash{variation} = $var_db_name;
      $species_hash{regulation} = $reg_db_name;
      
      # do we have SIFT or PolyPhen?
      if($var_db_name) {
        my $has_sift_poly = $self->has_sift_poly($dbc, $var_db_name, $species_id);
        $species_hash{$_} = $has_sift_poly->{$_} for keys %$has_sift_poly;
      }

      push @return, \%species_hash;
    }

    $species_count++;
    
    die("ERROR: Problem getting species and assembly names from $current_db_name; check coord_system table\n") unless $species_count;
  }
  
  return \@return;
}

sub is_strain {
  my ($self, $dbc, $current_db_name) = @_;

  my $sth = $dbc->prepare("select meta_value from ".$current_db_name.".meta where meta_key = 'species.strain';");
  $sth->execute();
  my $strain_value;
  $sth->bind_columns(\$strain_value);
  $sth->execute();
  $sth->fetch();
  $sth->finish();

  return 1 if $strain_value && $strain_value !~ /^reference/;

  return $current_db_name =~ /mus_musculus_.+?_(core|otherfeatures)/;
}

sub get_species_id_hash {
  my ($self, $dbc, $current_db_name) = @_;

  # get species names by id
  my $sth = $dbc->prepare("select species_id, meta_value from ".$current_db_name.".meta where meta_key = 'species.production_name';");
  $sth->execute();
  
  my ($species_id, $value, $species_ids);
  $sth->bind_columns(\$species_id, \$value);
  
  $species_ids->{$species_id} = $value while $sth->fetch();
  $sth->finish();

  return $species_ids;
}

sub get_assembly {
  my ($self, $dbc, $current_db_name, $species_id) = @_;

  my $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system WHERE species_id = ".$species_id." ORDER BY rank LIMIT 1;");
  $sth->execute();
  my $assembly;
  $sth->bind_columns(\$assembly);
  $sth->execute();
  $sth->fetch();
  $sth->finish();

  return $assembly;
}

sub has_var_db {
  my $self = shift;
  my $dbc = shift;
  my $current_db_name = shift;

  my $var_db_name = $current_db_name;
  $var_db_name =~ s/core|otherfeatures/variation/;
  my $has_var_db;
  my $sth = $dbc->prepare("SHOW DATABASES LIKE '$var_db_name';");
  $sth->execute();
  $sth->bind_columns(\$has_var_db);
  $sth->fetch;
  $sth->finish;

  return $has_var_db ? $var_db_name : undef;
}

sub has_sift_poly {
  my $self = shift;
  my $dbc = shift;
  my $var_db_name = shift;
  my $species_id = shift;
  $species_id ||= 0;

  my $sth = $dbc->prepare(qq{
    SELECT meta_key, meta_value
    FROM $var_db_name\.meta
    WHERE meta_key in ('sift_version','polyphen_version')
    AND (species_id = $species_id OR species_id IS NULL)
  });
  $sth->execute();
  
  my ($key, $val);
  $sth->bind_columns(\$key, \$val);
  
  my %data;

  while($sth->fetch) {
    $key =~ s/\_version//;
    $data{$key} = 'b';
  }
  $sth->finish();

  return \%data;
}

sub has_reg_build {
  my $self = shift;
  my $dbc = shift;
  my $current_db_name = shift;

  my $reg_db_name = $current_db_name;
  $reg_db_name =~ s/core|otherfeatures/funcgen/;
  my $has_reg_db;
  my $sth = $dbc->prepare("SHOW DATABASES LIKE '$reg_db_name';");
  $sth->execute();
  $sth->bind_columns(\$has_reg_db);
  $sth->fetch;
  $sth->finish;
  my $has_reg_build;

  if($has_reg_db) {
    $sth = $dbc->prepare("SELECT version FROM $reg_db_name.regulatory_build");
    $sth->execute();
    $sth->bind_columns(\$has_reg_build);
    $sth->fetch;
    $sth->finish;
  }

  return $has_reg_build ? $reg_db_name : undef;
}

1;

