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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::CreateDumpJobs;

use strict;
use warnings;

use File::Path qw(make_path);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub param_defaults {
  return {
    'merged' => 0,
    'eg' => 0,
    'core_jobs' => [],
    'refseq_jobs' => [],
    'variation_jobs' => [],
    'regulation_jobs' => [],
    'merge_jobs' => [],
  };
}

sub fetch_input {
  my $self = shift;

  my $group = $self->param('group');

  my %species_hash = map {$_ => $self->param($_)} qw(
    species
    assembly
    species_id
    is_multispecies
    variation
    regulation
    sift
    polyphen

    dbname
    host
    port
    user
    pass
    division
  );

  $species_hash{type} = $group eq 'otherfeatures' ? 'refseq' : 'core';

  # clearing the registry prevents a warning when we connect to
  # mutiple core DBs of the same species (e.g. core, otherfeatures)
  Bio::EnsEMBL::Registry->clear();

  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -group   => 'core',
    -species => $species_hash{species},
    -port    => $species_hash{port},
    -host    => $species_hash{host},
    -user    => $species_hash{user},
    -pass    => $species_hash{pass},
    -dbname  => $species_hash{dbname},
    -MULTISPECIES_DB => $species_hash{is_multispecies},
    -species_id => $species_hash{species_id},
  );

  if($species_hash{division}) {
    my $pipeline_dump_dir = $self->param('pipeline_dir')."/".$species_hash{division};
    $species_hash{pipeline_dump_dir} = $pipeline_dump_dir;
  }
  else {
    die "Can't find division for database ".$species_hash{dbname};
  }

  #Processing Collections
  my $mca = $dba->get_MetaContainerAdaptor;

  if($mca->is_multispecies == 1) {
    my $collection_db = $1 if($mca->dbc->dbname()=~/(.+)\_core/);
    $species_hash{dir_suffix} = "/".$collection_db;
    $species_hash{assembly}   = $mca->single_value_by_key('assembly.default');
  }
  $species_hash{db_version} = $mca->schema_version();

  # create join jobs
  my @join_jobs = (\%species_hash);

  if($group eq 'otherfeatures' && $self->param('merged')) {
    my %merged_job = %species_hash;
    $merged_job{type} = 'merged';
    push @join_jobs, \%merged_job;
    $self->param('merge_jobs', [\%merged_job]);
  }
  $self->param('join_jobs', \@join_jobs);

  # create per-chr jobs
  my $chr_jobs = $self->get_chr_jobs(\%species_hash, $dba);
  foreach my $type(qw(core refseq variation regulation)) {

    # get all the per-chr jobs of this type
    my @type_jobs = grep {$_->{type} eq $type} @$chr_jobs;    
    $self->param($type.'_jobs', \@type_jobs);
  }

  return;
}

sub write_output {
  my $self = shift;

  $self->dataflow_output_id($self->param('core_jobs'), 3);
  $self->dataflow_output_id($self->param('refseq_jobs'), 4);
  $self->dataflow_output_id($self->param('variation_jobs'), 5);
  $self->dataflow_output_id($self->param('regulation_jobs'), 6);

  # merge ens+refseq
  $self->dataflow_output_id($self->param('merge_jobs'), 7);

  # join, qc and finish can all use the same input_id
  $self->dataflow_output_id($self->param('join_jobs'), 8);
  $self->dataflow_output_id($self->param('join_jobs'), 9);
  $self->dataflow_output_id($self->param('join_jobs'), 10);
  
  return;
}

sub get_chr_jobs {
  my ($self, $species_hash, $dba) = @_;

  my $var_db = $self->param('group') eq 'core' ? $species_hash->{variation} : undef;
  my $reg_db = $self->param('group') eq 'core' ? $species_hash->{regulation} : undef;

  # get slices
  my $sa = $dba->get_SliceAdaptor;
  my @slices = @{$sa->fetch_all('toplevel')};
  push @slices, map {$_->alternate_slice} map {@{$_->get_all_AssemblyExceptionFeatures}} @slices;
  push @slices, @{$sa->fetch_all('lrg', undef, 1, undef, 1)} if $self->param('lrg') && $self->param('group') ne 'otherfeatures';

  # remove/sort out duplicates, in human you get 3 Y slices
  my %by_name;
  $by_name{$_->seq_region_name}++ for @slices;
  if(my @to_fix = grep {$by_name{$_} > 1} keys %by_name) {

    foreach my $name(@to_fix) {

      # remove those with duplicate name
      @slices = grep {$_->seq_region_name ne $name} @slices;

      # add a standard-fetched slice
      push @slices, $sa->fetch_by_region(undef, $name);
    }
  }

  # dump synonyms
  if($self->param('group') eq 'core') {
    make_path($species_hash->{pipeline_dump_dir}.'/synonyms');

    my $srsa = $dba->get_SeqRegionSynonymAdaptor();

    open SYN, sprintf(
      ">%s/synonyms/%s_%s_chr_synonyms.txt",
      $species_hash->{pipeline_dump_dir},
      $species_hash->{species},
      $species_hash->{assembly}
    ) or die "ERROR: Could not write to synonyms file\n $_";

    # synonyms can be indirect i.e. A <-> B <-> C
    # and there may not be a direct link between A <-> C in the DB
    # so let's allow for one level of indirection
    my $tree = {};
    
    my @all_syns = @{$srsa->fetch_all};

    # To prevent memory errors, species with many seq_regions will not search for indirect synonyms
    # With 4GB memory supplied, errors seem to start with species with ~50k seq regions, so 
    # we've chosen 40k as the limit for calculating indirect synonyms 
    my $syn_threshold = 40000; 
    foreach my $syn(@all_syns) {
      my $syn_slice = $sa->fetch_by_seq_region_id($syn->seq_region_id);
      next unless $syn_slice;
      my ($a, $b) = sort ($syn_slice->seq_region_name, $syn->name);
      $tree->{$a}->{$b} = 1;
      unless(scalar(@all_syns) > $syn_threshold){ 
        $tree->{$_}->{$b} = 1 for keys %{$tree->{$a} || {}};
        $tree->{$_}->{$a} = 1 for keys %{$tree->{$b} || {}};
      }
    }

    # now create uniq A <-> B / B <-> A pairs
    my %uniq;
    foreach my $a(keys %$tree) {
      foreach my $b(grep {$a ne $_} keys %{$tree->{$a}}) {
        $uniq{join("\t", sort ($a, $b))} = 1;
      }
    }

    print SYN "$_\n" for keys %uniq;

    close SYN;
  }

  # remove slices with no transcripts or variants on them
  # otherwise the dumper will create a load of "empty" cache files
  # in species with lots of unplaced scaffolds this means we create
  # masses of pointless directories which take ages to process
  my $ta = $dba->get_TranscriptAdaptor();

  @slices = grep {
    $ta->count_all_by_Slice($_) ||
    $self->count_vars($var_db, $species_hash, $_) ||
    $self->count_regfeats($reg_db, $species_hash, $_)
  } @slices;

  delete($self->{_var_dba}) if $self->{_var_dba};

  # now distribute the slices into jobs
  # jobs can contain multiple slices
  my @jobs;
  my $min_length = 10e6;
  my $added_length = 0;
  my %hash;

  foreach my $slice(sort {$b->end - $b->start <=> $a->end - $a->start} @slices) {
    unless(%hash) {
      %hash = %$species_hash;
      $hash{added_length} = 0;
    }

    push @{$hash{regions}}, {
      chr => $slice->seq_region_name,
      seq_region_id => $slice->get_seq_region_id,
      start => $slice->start,
      end => $slice->end,
    };

    $hash{added_length} += ($slice->end - $slice->start);

    if($hash{added_length} > $min_length) {
      push @jobs, @{$self->add_to_jobs(\%hash, $var_db, $reg_db)};
      $added_length = 0;
      %hash = ();
    }
  }

  push @jobs, @{$self->add_to_jobs(\%hash, $var_db, $reg_db)} if scalar keys %hash;

  # sort by type then length
  @jobs = sort {$a->{type} cmp $b->{type} || $b->{added_length} <=> $b->{added_length}} @jobs;

  return \@jobs;
}

sub count_vars {
  my ($self, $var_db_name, $species_hash, $slice) = @_;

  return 0 unless $var_db_name;

  unless(exists($self->{_var_dba})) {
    $self->{_var_dba} = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
      -group   => 'variation',
      -species => $species_hash->{species},
      -port    => $species_hash->{port},
      -host    => $species_hash->{host},
      -user    => $species_hash->{user},
      -pass    => $species_hash->{pass},
      -dbname  => $var_db_name,
      -MULTISPECIES_DB => $species_hash->{is_multispecies},
      -species_id => $species_hash->{species_id},
    );
  }

  my $sth = $self->{_var_dba}->dbc->prepare("SELECT COUNT(*) FROM variation_feature WHERE seq_region_id = ?");
  $sth->execute($slice->get_seq_region_id);

  my $v_count;
  $sth->bind_columns(\$v_count);
  $sth->fetch;
  $sth->finish;

  return $v_count;
}

sub count_regfeats {
  my ($self, $reg_db_name, $species_hash, $slice) = @_;

  return 0 unless $reg_db_name;

  unless(exists($self->{_rfa})) {
    my $dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
      -group   => 'funcgen',
      -species => $species_hash->{species},
      -port    => $species_hash->{port},
      -host    => $species_hash->{host},
      -user    => $species_hash->{user},
      -pass    => $species_hash->{pass},
      -dbname  => $reg_db_name,
      -MULTISPECIES_DB => $species_hash->{is_multispecies},
      -species_id => $species_hash->{species_id},
    );

    $self->{_rfa} = $dba->get_RegulatoryFeatureAdaptor;
  }

  return $self->{_rfa}->count_by_Slice_constraint($slice);
}

sub add_to_jobs {
  my ($self, $hash, $var_db, $reg_db) = @_;

  my @jobs;

  my %copy = %$hash;
  push @jobs, \%copy;

  if($var_db) {
    my %var = %{$hash};
    $var{type} = 'variation';
    push @jobs, \%var;
  }

  if($reg_db) {
    my %reg = %{$hash};
    $reg{type} = 'regulation';
    push @jobs, \%reg;
  }

  return \@jobs
}

1;

