=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy self the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, sselftware
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

# EnsEMBL module for Bio::EnsEMBL::VEP::Stats
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Stats - object for tracking VEP run stats

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Stats;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Variation::Utils::Constants;
use Bio::EnsEMBL::VEP::Constants;
use Bio::EnsEMBL::VEP::Utils qw(get_time);

# use CGI qw/:standard/;

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(no_stats)]);

  $self->{stats} = {counters => {}};

  return $self;
}

sub info {
  my $self = shift;
  $self->{info} = shift if @_;
  return $self->{info} ||= {};
}

sub start_time {
  my $self = shift;
  $self->{stats}->{run_time_start} = time();
  return $self->{stats}->{start_time} ||= get_time();
}

sub end_time {
  return $_[0]->{stats}->{end_time} ||= get_time();
}

sub run_time {
  return time() - $_[0]->{stats}->{run_time_start};
}

sub log_lines_read {
  $_[0]->{stats}->{lines_read} = $_[1];
}

sub log_fasta_chromosomes {
  my ($self, $fasta_db) = @_;

  return unless $fasta_db;
  
  $self->{stats}->{chr_lengths} ||= {
    map {$_ => $fasta_db->length($_)}
    $fasta_db->isa('Bio::DB::Fasta') ?
    $fasta_db->get_all_primary_ids :
    $fasta_db->get_all_sequence_ids
  };
}

sub log_db_chromosomes {
  $_[0]->{stats}->{chr_lengths} ||= $_[1];
}

sub log_VariationFeature {
  my $self = shift;
  my $vf = shift;
  my $hash = shift;

  return if $self->{no_stats};

  my $stats = $self->{stats}->{counters};

  # basic count
  $stats->{var_count}++;

  # position
  $stats->{chr}->{$vf->{chr}}->{1e6 * int($vf->{start} / 1e6)}++;

  # most severe consequence
  $stats->{var_cons}->{$vf->display_consequence}++;

  # known variants
  $stats->{existing}++ if $vf->{existing};

  # get stats
  my $so_term = $vf->class_SO_term;
  if(defined($so_term)) {
    $stats->{classes}->{$so_term}++;
    $stats->{allele_changes}->{$vf->allele_string}++ if $so_term eq 'SNV';
  }
}

sub log_VariationFeatureOverlapAllele {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  return if $self->{no_stats};

  my $stats = $self->{stats}->{counters};

  # consequence type(s)
  $stats->{consequences}->{$_}++ for @{$hash->{Consequence}};

  # feature types
  $stats->{lc($hash->{Feature_type})}->{$hash->{Feature}}++ if $hash->{Feature};

  # run sub-type method if available
  my $method = 'log_'.(split('::', ref($vfoa)))[-1];
  $self->$method(@_) if $self->can($method);
}

sub log_TranscriptVariationAllele {
  my $self = shift;
  $self->log_BaseTranscriptVariationAllele(@_);
}

sub log_TranscriptStructuralVariationAllele {
  my $self = shift;
  $self->log_BaseTranscriptVariationAllele(@_);
}

sub log_BaseTranscriptVariationAllele {
  my $self = shift;
  my ($vfoa, $hash) = @_;

  my $stats = $self->{stats}->{counters};

  # protein pos
  if($vfoa->_pre_consequence_predicates->{coding}) {
    my $tr = $vfoa->feature;
    my $protein_length =
      $tr->{_variation_effect_feature_cache}->{peptide} ?
      length($tr->{_variation_effect_feature_cache}->{peptide}) :
      $tr->translation->length;

    if(my $protein_pos = $vfoa->base_transcript_variation->translation_start) {
      $stats->{protein_pos}->{int(10 * ($protein_pos / $protein_length))}++ if $protein_length;
    }
  }

  # gene counts
  $stats->{gene}->{$hash->{Gene}}++;
}

sub log_sift_polyphen {
  my ($self, $tool, $pred) = @_;
  $self->{stats}->{counters}->{$tool}->{$pred}++;
}

sub finished_stats {
  my $self = shift;

  if(!exists($self->{finished_stats})) {

    my $stats = $self->{stats};

    # merge in counters data
    my $counters = $stats->{counters};
    $stats->{$_} = $counters->{$_} for keys %$counters;

    # convert gene and transcript hashes to counts
    for my $type(qw(gene transcript regulatoryfeature)) {
      $stats->{$type} = scalar keys %{$stats->{$type}} if defined $stats->{$type};
    }
    
    # tot up chromosome counts
    foreach my $chr(keys %{$stats->{chr}}) {
      $stats->{chr_totals}->{$chr} += $stats->{chr}->{$chr}->{$_} for keys %{$stats->{chr}->{$chr}};
      
      my $start = 0;
      my %tmp;
      
      while($start <= $stats->{chr_lengths}->{$chr}) {
        $tmp{$start / 1e6} = $stats->{chr}->{$chr}->{$start} || 0;
        $start += 1e6;
      }
      
      $stats->{chr}->{$chr} = \%tmp;
    }
    
    # convert allele changes to Ts/Tv
    my $ts_tv = \%Bio::EnsEMBL::VEP::Constants::TS_TV;
    map {$stats->{ts_tv}->{$ts_tv->{$_}} += $stats->{allele_changes}->{$_}} keys %{$stats->{allele_changes}} if $stats->{allele_changes};
    
    # flesh out protein_pos
    if(defined($stats->{protein_pos})) {
      if(defined($stats->{protein_pos}->{10})) {
        $stats->{protein_pos}->{9} += $stats->{protein_pos}->{10};
        delete $stats->{protein_pos}->{10};
      }
      $stats->{protein_pos}->{$_} ||= 0 for (0..9);
      
      my %tmp = map {$_.'0-'.($_+1).'0%' => $stats->{protein_pos}->{$_}} keys %{$stats->{protein_pos}};
      $stats->{protein_pos} = \%tmp;
    }
    
    # coding cons
    foreach my $con(
      map {$_->SO_term}
      grep {$_->{include} && $_->{include}->{coding}}
      values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES
    ) {
      $stats->{coding}->{$con} = $stats->{consequences}->{$con} if $stats->{consequences}->{$con};
    }

    $self->{finished_stats} = {
      charts => $self->generate_chart_data($stats),
      run_stats => $self->generate_run_stats($stats),
      general_stats => $self->generate_general_stats($stats),
    }
  }

  return $self->{finished_stats};
}

sub generate_chart_data {
  my $self = shift;
  my $stats = shift;

  # get some data from Constants
  my $colour_keys = \%Bio::EnsEMBL::VEP::Constants::COLOUR_KEYS;
  my %cons_ranks =
    map { $_->{SO_term} => $_->{rank} }
    values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;

  # create pie chart hashes
  my @charts = (
    {
      id => 'var_class',
      title => 'Variant classes',
      header => ['Variant class', 'Count'],
      data => $stats->{classes},
      type => 'pie',
      sort => 'value',
      height => 200,
    },
    {
      id => 'var_cons',
      title => 'Consequences (most severe)',
      header => ['Consequence type', 'Count'],
      data => $stats->{var_cons},
      type => 'pie',
      sort => \%cons_ranks,
      colours => $colour_keys->{consequences},
    },
    {
      id => 'consequences',
      title => 'Consequences (all)',
      header => ['Consequence type', 'Count'],
      data => $stats->{consequences},
      type => 'pie',
      sort => \%cons_ranks,
      colours => $colour_keys->{consequences},
    },
    {
      id => 'coding',
      title => 'Coding consequences',
      header => ['Consequence type', 'Count'],
      data => $stats->{coding},
      type => 'pie',
      sort => \%cons_ranks,
      colours => $colour_keys->{consequences},
    }
  );
  
  foreach my $tool(qw(SIFT PolyPhen)) {
    my $lc_tool = lc($tool);
    
    push @charts, {
      id => $lc_tool,
      title => $tool.' summary',
      header => ['Prediction', 'Count'],
      data => $stats->{$tool},
      type => 'pie',
      height => 200,
      sort => 'value',
      colours => $colour_keys->{$lc_tool},
    } if $self->param($lc_tool);
  }
  
  push @charts, {
    id => 'chr',
    title => 'Variants by chromosome',
    header => ['Chromosome','Count'],
    data => $stats->{chr_totals},
    sort => 'chr',
    type => 'bar',
    options => '{legend: {position: "none"}}',
  };
  
  foreach my $chr(sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %{$stats->{chr}}) {
    my $chr_id = $chr;
    $chr_id =~ s/\./\_/g;
    
    push @charts, {
      id => 'chr_'.$chr_id,
      title => 'Distribution of variants on chromosome '.$chr,
      header => ['Position (mb)', 'Count'],
      data => $stats->{chr}->{$chr},
      sort => 'chr',
      type => 'area',
      options => '{hAxis: {title: "Position (mb)", textStyle: {fontSize: 8}}, legend: {position: "none"}}',
      no_table => 1,
      no_link => 1,
    };
  }
  
  push @charts, {
    id => 'protein',
    title => 'Position in protein',
    header => ['Position in protein (percentile)','Count'],
    data => $stats->{protein_pos},
    sort => 'chr',
    type => 'bar',
    no_table => 1,
    options => '{hAxis: {title: "Position in protein (percentile)", textStyle: {fontSize: 10}}, legend: {position: "none"}}',
  };

  return \@charts;
}

sub generate_run_stats {
  my $self = shift;
  my $stats = shift;

  my $info = $self->{info};
    
  return [
    [sprintf('VEP version (API) %i (%i)', $info->{vep_version}, $info->{api_version})],
    ['Cache/Database', ($info->{cache_dir} ? $info->{cache_dir} : sprintf('%s on %s', $info->{db_name}, $info->{db_host}))],
    ['Species', $self->species],
    ['Command line options', '<pre>'.join(" ", %{$self->config->_raw_config}).'</pre>'],
    ['Start time', $self->start_time],
    ['End time', $self->end_time],
    ['Run time', $self->run_time." seconds"],
    ['Input file', $self->param('input_file')],
    [
      'Output file',
      $self->param('output_file')#.
      # (defined($config->{html}) ? ' '.a({href => $config->{output_file}.'.html'}, '[HTML]') : '').
      # ' '.a({href => $config->{output_file}}, '[text]')
    ],
  ];
}

sub generate_general_stats {
  my $self = shift;
  my $stats = shift;

  return [
    ['Lines of input read', $stats->{lines_read}],
    ['Variants processed', $stats->{var_count}],
    # ['Variants remaining after filtering', $stats->{filter_count}],
    # ['Lines of output written', $stats->{out_count}],
    [
      'Novel / existing variants',
      defined($stats->{existing}) ?
      sprintf("%s (%.1f\%) / %s (%.1f\%)",
        $stats->{var_count} - $stats->{existing},
        100 * (($stats->{var_count} - $stats->{existing}) / $stats->{var_count}),
        $stats->{existing},
        100 * ($stats->{existing} / $stats->{var_count}),
      )
      : '-'
    ],
    ['Overlapped genes', $stats->{gene}],
    ['Overlapped transcripts', $stats->{transcript}],
    ['Overlapped regulatory features', $stats->{regulatoryfeature} || '-'],
  ];
}

sub sort_keys {
  my $self = shift;
  my $data = shift;
  my $sort = shift;
  
  my @keys;
  
  # sort data
  if(defined($sort)) {
    if($sort eq 'chr') {
      @keys = sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %{$data};
    }
    elsif($sort eq 'value') {
      @keys = sort {$data->{$a} <=> $data->{$b}} keys %{$data};
    }
    elsif(ref($sort) eq 'HASH') {
      @keys = sort {$sort->{$a} <=> $sort->{$b}} keys %{$data};
    }
  }
  else {
    @keys = keys %{$data};
  }
  
  return \@keys;
}

sub dump_text {
  my $self = shift;
  my $fh = shift || *STDERR;

  my $finished_stats = $self->finished_stats;

  print $fh "[VEP run statistics]\n";
  print $fh join("\t", map {s/\<.+?\>//g; $_} @{$_})."\n" for @{$finished_stats->{run_stats}};
  
  print $fh "\n[General statistics]\n";
  print $fh join("\t", map {s/\<.+?\>//g; $_} @{$_})."\n" for @{$finished_stats->{general_stats}};
  
  foreach my $chart(@{$finished_stats->{charts}}) {
    print $fh "\n[".$chart->{title}."]\n";
    print $fh join("\t", ($_, $chart->{data}->{$_}))."\n" for @{$self->sort_keys($chart->{data}, $chart->{sort})};
  }
}

sub dump_html {
  my $self = shift;
  my $fh = shift;
}

1;