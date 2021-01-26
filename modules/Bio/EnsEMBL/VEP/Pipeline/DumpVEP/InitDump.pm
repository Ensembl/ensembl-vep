=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
use File::Path qw(mkpath);

my $DEBUG = 0;

sub param_defaults {
  return {
    'group' => 'core'
  };
}

sub run {
  my $self = shift;
  my $species = $self->required_param('species');
  my $group = $self->required_param('group');

  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
  my $dbc = $dba->dbc();
  my $current_db_name = $dbc->dbname();
 
  my $var_db_name = $self->has_var_db($dbc, $current_db_name);
  my $dbc_var;

  if($var_db_name){
    my $dba_var = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'variation');
    $dbc_var = $dba_var->dbc();
  }

  #Special case for otherfeatures
  if ($current_db_name =~ /otherfeatures/) {
    if ($self->has_refseq($dbc, $current_db_name) <= 0) {
      $self->warning("$current_db_name doesn't have any RefSeq Transcripts, skipping it");
      return;
    }
  }
  
  mkpath($self->param('pipeline_dir')) unless -d $self->param('pipeline_dir');

  $self->param(
    'species_jobs',
    [
      map {@{$self->generate_species_jobs($_,$group,$dba,$dbc,$dbc_var,$current_db_name)}}
      $species,
    ]
  );

  $self->dataflow_output_id($self->param('species_jobs'), 2);
  return;
}

sub generate_species_jobs {
  my ($self, $species, $group, $dba, $dbc, $dbc_var, $current_db_name) = @_;

  my @return;
  
  my $version = $self->param('eg_version') || $self->required_param('ensembl_release');

  my $species_id = $self->get_species_id($dbc, $current_db_name, $species);
    
  # do we have a variation DB?
  my $var_db_name = $self->has_var_db($dbc, $current_db_name);
    
  # do we have a regulation DB?
  my $reg_db_name = $self->has_reg_build($dbc, $current_db_name);

  my $species_count = 0;
    
  my $assembly = $self->get_assembly($dbc, $current_db_name, $species_id);
  next unless $assembly;
      
  # copy server details
  my %species_hash;
      
  $species_hash{species} = $species;
  $species_hash{species_id} = $species_id;
  $species_hash{assembly} = $assembly;
  $species_hash{dbname} = $current_db_name;
  $species_hash{group} = $group;
  $species_hash{is_multispecies} = $current_db_name =~ "_collection" ? 1 : 0;
  $species_hash{variation} = $var_db_name;
  $species_hash{regulation} = $reg_db_name;
  $species_hash{host} = $dbc->host();
  $species_hash{port} = $dbc->port();
  $species_hash{user} = $dbc->username();
  $species_hash{pass} = $dbc->password();
  $species_hash{division} = $self->division($dba);

  # do we have SIFT or PolyPhen?
  if($var_db_name) {
    my $has_sift_poly = $self->has_sift_poly($dbc, $var_db_name, $species_id);
    $species_hash{$_} = $has_sift_poly->{$_} for keys %$has_sift_poly;
    
    $self->pubmed($dbc_var, $current_db_name, $species_id);
    $self->synonyms($dbc_var, $current_db_name, $species_id);
    
  }

  push @return, \%species_hash;

  $species_count++;

  die("ERROR: Problem getting species and assembly names from $current_db_name; check coord_system table\n") unless $species_count;

  $dbc->disconnect_if_idle();

  return \@return;
}

sub pubmed {
  my ($self, $dbc, $current_db_name, $species_id) = @_;
  $species_id ||= 0;

  my %pm;
  my $pipeline_dump_dir = $self->param('pipeline_dir');

  my $file = sprintf(
      '%s/pubmed_%s_%s_%s.txt',
      $pipeline_dump_dir,
      $self->required_param('species'),
      $self->required_param('ensembl_release'),
      $self->get_assembly($dbc, $current_db_name, $species_id)
      );
      
  my $lock = $file.'.lock';
  unlink($file) if -e $file;
  unlink($lock) if -e $lock;
  my $sleep_count = 0;
  if(-e $lock) {
    while(-e $lock) {
      sleep 1;
      die("I've been waiting for $lock to be removed for $sleep_count seconds, something may have gone wrong\n") if ++$sleep_count > 900;
    }
  }
  if(-e $file) {
    unlink($file) or die "Failed to delete pubmed file $file $!";
  }
    
  if($dbc) {
    open OUT, ">$lock" or die "Could not write to pubmed lock file $lock $!";
    print OUT "$$\n";
    close OUT;
    $self->{_lock_file} = $lock;

    my $sth = $dbc->prepare(qq{
      SELECT v.name, GROUP_CONCAT(p.pmid)
      FROM variation v, variation_citation c, publication p
      WHERE v.variation_id = c.variation_id
      AND c.publication_id = p.publication_id
      AND p.pmid IS NOT NULL
      GROUP BY v.variation_id
    });
    $sth->execute();

    my ($v, $p);
    $sth->bind_columns(\$v, \$p);

    open OUT, ">$file" or die "Could not write to pubmed file $file $!";
    while($sth->fetch()) {
      $pm{$v} = $p;
      print OUT "$v\t$p\n";
    }

    close OUT;

    unlink($lock) or die "Failed to delete pubmed lock file $lock $!";

    $sth->finish();
  } 
  return \%pm;
}

sub synonyms {
  my ($self, $dbc, $current_db_name, $species_id) = @_;
  $species_id ||= 0;

  my %pm;
  my $pipeline_dump_dir = $self->param('pipeline_dir');

  my $file = sprintf(
      '%s/var_synonyms_%s_%s_%s.txt',
      $pipeline_dump_dir,
      $self->required_param('species'),
      $self->required_param('ensembl_release'),
      $self->get_assembly($dbc, $current_db_name, $species_id)
      );

  my $lock = $file.'.lock';
  unlink($file) if -e $file;
  unlink($lock) if -e $lock;
  my $sleep_count = 0;
  if(-e $lock) {
    while(-e $lock) {
      sleep 1;
      die("I've been waiting for $lock to be removed for $sleep_count seconds, something may have gone wrong\n") if ++$sleep_count > 900;
    }
  }
  if(-e $file) {
    unlink($file) or die "Failed to delete var synonyms file $file $!";
  }

  if($dbc) {
    open OUT, ">$lock" or die "Could not write to synonyms lock file $lock $!";
    print OUT "$$\n";
    close OUT;
    $self->{_lock_file} = $lock;

    my $sth = $dbc->prepare(qq{
    select subtab.variation_id, group_concat(subtab.str separator '--') as str 
    from (
	   select variation_id, concat(s.name, '::', group_concat(vs.name order by vs.variation_id separator ',')) as str 
	   from variation_synonym vs 
	   inner join source s on s.source_id = vs.source_id 
	   where s.source_id in ( 
               select source_id 
	       from source 
	       where name not like '%dbSNP%')
	    group by variation_id, s.name) as subtab 
    group by subtab.variation_id;
    });
    $sth->execute();

    my ($v, $p);
    $sth->bind_columns(\$v, \$p);

    open OUT, ">$file" or die "Could not write to synonyms file $file $!";
    while($sth->fetch()) {
      $pm{$v} = $p;
      print OUT "$v\t$p\n";
    }

    close OUT;

    unlink($lock) or die "Failed to delete var synonyms lock file $lock $!";

    $sth->finish();
  }
  return \%pm;
}

sub get_species_id {
  my ($self, $dbc, $current_db_name, $species) = @_;

  # get species id
  my $species_id = $dbc->sql_helper()
                    ->execute_simple( -SQL =>qq/select species_id from $current_db_name.meta where meta_key = 'species.production_name' and meta_value ='$species';/);
 
  return $species_id->[0];
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
  my ($self, $dbc, $current_db_name) = @_;

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

sub has_refseq {
  my ($self, $dbc, $current_db_name) = @_;

  my $sth = $dbc->prepare(qq{
      SELECT COUNT(*)
      FROM $current_db_name\.transcript t LEFT JOIN $current_db_name\.xref x ON x.xref_id = t.display_xref_id
      WHERE t.source LIKE '%refseq%'
      OR x.display_label like 'NM%'
  });
  
  $sth->execute;
    
  my $count;
  $sth->bind_columns(\$count);
  $sth->fetch;
  $sth->finish();

  return $count ? $count : undef;
}

sub has_sift_poly {
  my ($self, $dbc, $var_db_name, $species_id) = @_;
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
  my ($self, $dbc, $current_db_name) = @_;

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

sub division {
    my ($self, $dba) = @_;
    my ($division) = @{$dba->get_MetaContainer()->list_value_by_key('species.division')};
    return if ! $division;
    $division =~ s/^Ensembl//;

return lc($division);
}

1;

