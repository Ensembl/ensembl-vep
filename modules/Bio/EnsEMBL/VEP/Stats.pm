=head1 LICENSE

Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

=head1 SYNOPSIS

my $stats = Bio::EnsEMBL::VEP::Haplo::Stats->new({
  config => $config
});

my $start_time = $stats->start_time();

$stats->log_VariationFeature($vf);

my $end_time = $stats->end_time();

=head1 DESCRIPTION

The Stats class is used to track all runtime and data statistics
for a VEP run.

It also generates summary statistics output, by default in an HTML
report file containing tables and charts.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::Stats;

use base qw(Bio::EnsEMBL::VEP::BaseVEP);

use Bio::EnsEMBL::Variation::Utils::Constants;
use Bio::EnsEMBL::VEP::Constants;
use Bio::EnsEMBL::VEP::Utils qw(get_time);


=head2 new

  Arg 1      : hashref $args
               {
                 config => Bio::EnsEMBL::VEP::Config,
               }
  Example    : $stats = Bio::EnsEMBL::VEP::Stats->new({config => $config});
  Description: Create a new Bio::EnsEMBL::VEP::Stats object.
  Returntype : Bio::EnsEMBL::VEP::Stats
  Exceptions : none
  Caller     : BaseVEP
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my $self = $class->SUPER::new(@_);

  # add shortcuts to these params
  $self->add_shortcuts([qw(no_stats)]);

  $self->{stats} = {counters => {}};

  return $self;
}


=head2 info

  Arg 1      : (optional) hashref $info
  Example    : $info = $stats->info();
  Description: Get/set the info hashref. Used to store version and annotation
               source data, as retrieved by get_output_header_info() in Runner.
  Returntype : hashref
  Exceptions : none
  Caller     : Runner::init()
  Status     : Stable

=cut

sub info {
  my $self = shift;
  $self->{info} = shift if @_;
  return $self->{info} ||= {};
}


=head2 start_time

  Example    : $time = $stats->start_time();
  Description: Get start timestamp of this VEP run. Set to current time the first
               time this method is called.
  Returntype : string
  Exceptions : none
  Caller     : Runner::init()
  Status     : Stable

=cut

sub start_time {
  my $self = shift;
  $self->{stats}->{run_time_start} ||= time();
  return $self->{stats}->{start_time} ||= get_time();
}


=head2 end_time

  Example    : $time = $stats->end_time();
  Description: Get end timestamp of this VEP run. Set to current time the first
               time this method is called.
  Returntype : string
  Exceptions : none
  Caller     : Runner::dump_stats()
  Status     : Stable

=cut

sub end_time {
  return $_[0]->{stats}->{end_time} ||= get_time();
}


=head2 run_time

  Example    : $runtime = $stats->run_time();
  Description: Get run time in seconds (current time - start time)
  Returntype : int
  Exceptions : none
  Caller     : generate_run_stats()
  Status     : Stable

=cut

sub run_time {
  return time() - $_[0]->{stats}->{run_time_start};
}


=head2 log_lines_read

  Arg 1      : int $num_lines
  Example    : $stats->log_lines_read($num_lines);
  Description: Store the total number of lines read, used in stats output
  Returntype : none
  Exceptions : none
  Caller     : InputBuffer
  Status     : Stable

=cut

sub log_lines_read {
  $_[0]->{stats}->{lines_read} = $_[1];
}


=head2 log_chromosomes

  Arg 1      : hashref $chr_names_and_lengths
  Example    : $stats->log_chromosomes($lengths);
  Description: Logs the chromosome names and lengths as retrieved from
               core database or FASTA db
  Returntype : none
  Exceptions : none
  Caller     : BaseVEP
  Status     : Stable

=cut

sub log_chromosomes {
  $_[0]->{stats}->{chr_lengths} ||= $_[1];
}


=head2 log_VariationFeature

  Arg 1      : Bio::EnsEMBL::Variation::BaseVariationFeature $vf
  Arg 2      : hashref $vf_hash
  Example    : $stats->log_VariationFeature($vf);
  Description: Log statistics for a VariationFeature object
  Returntype : none
  Exceptions : none
  Caller     : OutputFactory
  Status     : Stable

=cut

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

    if($so_term eq 'SNV') {
      my @alleles = split('/', $vf->allele_string);
      my $ref_allele = shift @alleles;
      $stats->{allele_changes}->{$ref_allele.'/'.$_}++ for @alleles;
    }
  }
}


=head2 log_VariationFeatureOverlapAllele

  Arg 1      : Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele $vfoa
  Arg 2      : hashref $vf_hash
  Example    : $stats->log_VariationFeatureOverlapAllele($vfoa);
  Description: Log statistics for a VariationFeatureOverlapAllele object
  Returntype : none
  Exceptions : none
  Caller     : OutputFactory
  Status     : Stable

=cut

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


=head2 log_TranscriptVariationAllele

  Arg 1      : Bio::EnsEMBL::Variation::TranscriptVariationAllele $tva
  Arg 2      : hashref $vf_hash
  Example    : $stats->log_TranscriptVariationAllele($tva);
  Description: Log statistics for a TranscriptVariationAllele object
  Returntype : none
  Exceptions : none
  Caller     : OutputFactory
  Status     : Stable

=cut

sub log_TranscriptVariationAllele {
  my $self = shift;
  $self->log_BaseTranscriptVariationAllele(@_);
}


=head2 log_TranscriptStructuralVariationAllele

  Arg 1      : Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele $tsva
  Arg 2      : hashref $vf_hash
  Example    : $stats->log_TranscriptStructuralVariationAllele($tsva);
  Description: Log statistics for a TranscriptStructuralVariationAllele object
  Returntype : none
  Exceptions : none
  Caller     : OutputFactory
  Status     : Stable

=cut

sub log_TranscriptStructuralVariationAllele {
  my $self = shift;
  $self->log_BaseTranscriptVariationAllele(@_);
}


=head2 log_BaseTranscriptVariationAllele

  Arg 1      : Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele $btva
  Arg 2      : hashref $vf_hash
  Example    : $stats->log_BaseTranscriptVariationAllele($btva);
  Description: Log statistics for a BaseTranscriptVariationAllele object
  Returntype : none
  Exceptions : none
  Caller     : log_TranscriptVariationAllele()
               log_TranscriptStructuralVariationAllele()
  Status     : Stable

=cut

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


=head2 log_sift_polyphen

  Arg 1      : string $tool_name
  Arg 2      : string $prediction
  Example    : $stats->log_sift_polyphen('SIFT', 'deleterious');
  Description: Log statistics for a SIFT or PolyPhen prediction
  Returntype : none
  Exceptions : none
  Caller     : OutputFactory
  Status     : Stable

=cut

sub log_sift_polyphen {
  my ($self, $tool, $pred) = @_;
  $self->{stats}->{counters}->{$tool}->{$pred}++;
}


=head2 increment_filtered_variants

  Arg 1      : int $count
  Example    : $stats->increment_filtered_variants(3);
  Description: Log a number of filtered variants
  Returntype : none
  Exceptions : none
  Caller     : AnnotationSource::Cache::BaseCacheVariation
  Status     : Stable

=cut

sub increment_filtered_variants {
  my ($self, $count) = @_;
  $self->{stats}->{counters}->{filtered_variants} += $count;
}


=head2 finished_stats

  Example    : $finished = $stats->finished_stats();
  Description: Finalise and return finished stats hashref
  Returntype : hashref
  Exceptions : none
  Caller     : dump_text(), dump_html()
  Status     : Stable

=cut

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
    if($stats->{chr_lengths}) {
      foreach my $chr(keys %{$stats->{chr}}) {
        $stats->{chr_totals}->{$chr} += $stats->{chr}->{$chr}->{$_} for keys %{$stats->{chr}->{$chr}};
        
        unless($stats->{chr_lengths}->{$chr}) {
          delete $stats->{chr}->{$chr};
          next;
        }

        my $start = 0;
        my %tmp;
        
        while($start <= $stats->{chr_lengths}->{$chr}) {
          $tmp{$start / 1e6} = $stats->{chr}->{$chr}->{$start} || 0;
          $start += 1e6;
        }
        
        $stats->{chr}->{$chr} = \%tmp;
      }
    }
    
    # convert allele changes to Ts/Tv
    my $ts_tv = \%Bio::EnsEMBL::VEP::Constants::TS_TV;
    map {$stats->{ts_tv}->{$ts_tv->{$_}} += $stats->{allele_changes}->{$_}} grep {$ts_tv->{$_}} keys %{$stats->{allele_changes}} if $stats->{allele_changes};
    
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


=head2 generate_chart_data

  Arg 1      : hashref $stats
  Example    : $chart_data = $stats->generate_chart_data($stats);
  Description: Generates chart configurations and data based on available
               stats.
  Returntype : arrayref
  Exceptions : none
  Caller     : finished_stats()
  Status     : Stable

=cut

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
  } if $stats->{chr_totals};
  
  foreach my $chr(sort {($a !~ /^\d+$/ || $b !~ /^\d+/ || $a =~ /^\d\w/ || $b =~ /^\d\w/ ) ? $a cmp $b : $a <=> $b} keys %{$stats->{chr}}) {
    my $chr_id = $chr;
    $chr_id =~ s/\.|-/\_/g;

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
  } if $stats->{protein_pos};

  return \@charts;
}


=head2 generate_run_stats

  Arg 1      : hashref $stats
  Example    : $run_stats = $stats->generate_run_stats($stats);
  Description: Generates VEP run stats.
  Returntype : arrayref
  Exceptions : none
  Caller     : finished_stats()
  Status     : Stable

=cut

sub generate_run_stats {
  my $self = shift;
  my $stats = shift;

  my $info = $self->{info};

  # process command line opts
  my %raw = %{$self->config->_raw_config};
  my @opts;
  foreach my $key(sort keys %raw) {
    my $val = $raw{$key};
    if(ref($val) eq 'ARRAY') {
      next unless scalar @$val;
      push @opts, '--'.$key, $_ for @$val;
    }
    else {
      next if $val eq '0';
      push @opts, '--'.$key;
      push @opts, $val unless $val eq '1';
    }
  }

  my @return = (
    ['Species', $self->species],
    ['Command line options', '<pre>'.join(" ", @opts).'</pre>'],
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
  );

  my @cache_db_strings;
  if($info->{cache_dir}) {
    push @cache_db_strings, "Cache: ".$info->{cache_dir};
  }
  if($self->param('database') or ($self->param('cache') && !$self->param('offline'))) {
    push @cache_db_strings, sprintf('%s on %s', $info->{db_name}, $info->{db_host});
  }
  foreach my $custom(@{$info->{custom_info} || []}) {
    push @cache_db_strings, sprintf('Custom: %s (%s)', $custom->{file}, $custom->{type});
  }

  unshift @return, ['Annotation sources', join("; ", @cache_db_strings)];
  unshift @return, ['VEP version (API)', sprintf(' %i (%i)', $info->{vep_version}, $info->{api_version})];

  return \@return;
}


=head2 generate_general_stats

  Arg 1      : hashref $stats
  Example    : $run_stats = $stats->generate_run_stats($stats);
  Description: Generates general VEP stats, not run-related but not
               suitable for charting e.g. no. overlapped genes.
  Returntype : arrayref
  Exceptions : none
  Caller     : finished_stats()
  Status     : Stable

=cut

sub generate_general_stats {
  my $self = shift;
  my $stats = shift;

  return [
    ['Lines of input read', $stats->{lines_read}],
    ['Variants processed', $stats->{var_count}],
    ['Variants filtered out', $stats->{filtered_variants} || 0],
    # ['Lines of output written', $stats->{out_count}],
    [
      'Novel / existing variants',
      defined($stats->{existing}) ?
      sprintf('%s (%.1f) / %s (%.1f)',
        $stats->{var_count} - $stats->{existing},
        100 * (($stats->{var_count} - $stats->{existing}) / $stats->{var_count}),
        $stats->{existing},
        100 * ($stats->{existing} / $stats->{var_count}),
      )
      : '-'
    ],
    ['Overlapped genes', $stats->{gene} || 0],
    ['Overlapped transcripts', $stats->{transcript} || 0],
    ['Overlapped regulatory features', $stats->{regulatoryfeature} || ($self->param('regulatory') ? 0 : '-')],
  ];
}


=head2 sort_keys

  Arg 1      : hashref $data
  Arg 2      : scalar $sort_type
  Example    : $run_stats = $stats->generate_run_stats($stats);
  Description: Sorts the keys from a hashref of data according to various strategies
               as set by $sort_type:
                - string 'chr'   : sort alphanumerically
                - string 'value' : sort numerically by values in $data
                - hashref        : sort numerically by values of keys in $sort
  Returntype : arrayref of strings
  Exceptions : none
  Caller     : dump_text(), stats_html_head()
  Status     : Stable

=cut

sub sort_keys {
  my $self = shift;
  my $data = shift;
  my $sort = shift;
  
  my @keys;
  
  # sort data
  if(defined($sort)) {
    if($sort eq 'chr') {
      @keys = sort {($a !~ /^\d+$/ || $b !~ /^\d+/ || $a =~ /^\d\w/ || $b =~ /^\d\w/ ) ? $a cmp $b : $a <=> $b} keys %{$data};
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


=head2 dump_text

  Arg 1      : (optional) filehandle $fh
  Example    : $stats->dump_text(*STDOUT);
  Description: Dumps statistics as text to filehandle; STDERR is used if no
               filehandle given.
  Returntype : bool
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub dump_text {
  my $self = shift;
  my $fh = shift || *STDERR;

  my $finished_stats = $self->finished_stats;

  print $fh "[VEP run statistics]\n";
  print $fh join("\t", map {s/\<.+?\>//g; $_} @{$_})."\n" for @{$finished_stats->{run_stats}};
  
  print $fh "\n[General statistics]\n";
  print $fh join("\t", map {s/\<.+?\>//g; $_} grep {defined($_)} @{$_})."\n" for @{$finished_stats->{general_stats}};
  
  foreach my $chart(@{$finished_stats->{charts}}) {
    print $fh "\n[".$chart->{title}."]\n";
    print $fh join("\t", ($_, $chart->{data}->{$_}))."\n" for @{$self->sort_keys($chart->{data}, $chart->{sort})};
  }

  return 1;
}


=head2 dump_html

  Arg 1      : filehandle $fh
  Example    : $stats->dump_html($fh);
  Description: Dumps statistics as HTML to filehandle
  Returntype : bool
  Exceptions : none
  Caller     : Runner
  Status     : Stable

=cut

sub dump_html {
  my $self = shift;
  my $fh = shift;

  my $finished_stats = $self->finished_stats;

  print $fh $self->stats_html_head($finished_stats->{charts});
      
  # create menu
  print $fh
    '<div class="sidemenu">'.
      '<div class="sidemenu_head">Links</div>'.
      '<div class="sidemenu_body">'.
        '<ul>'.
          join('', map {sprintf('<li><a href="#%s">%s</a></li>', $_->[0], $_->[1])} (
            ['masthead', 'Top of page'],
            ['run_stats', 'VEP run statistics'],
            ['gen_stats', 'General statistics'],
            map {
              [$_->{id}, $_->{title}]
            } grep { !$_->{no_link} } @{$finished_stats->{charts}},
          )).
        '</ul>'.
      '</div>'.
    '</div>';
  
  print $fh "<div class='main_content'>";
  
  print $fh
    '<h3 id="run_stats">VEP run statistics</h3>'.
    '<table class="stats_table">'.
      join('', map {'<tr>'.join('', map {'<td>'.$_.'</td>'} @$_).'</tr>'} @{$finished_stats->{run_stats}}).
    '</table>';
  
  # vars in/out stats
  print $fh
    '<h3 id="gen_stats">General statistics</h3>'.
    '<table class="stats_table">'.
      join('', map {'<tr>'.join('', map {'<td>'.($_ || 0).'</td>'} @$_).'</tr>'} @{$finished_stats->{general_stats}}).
    '</table>';
  
  foreach my $chart(@{$finished_stats->{charts}}) {
    my $height = $chart->{height} || ($chart->{type} eq 'pie' ? '400' : '200');
    
    print $fh
      '<hr/>'.
      sprintf('<h3 id="%s">%s</h3>', $chart->{id}, $chart->{title}).
      sprintf('<div id="%s" style="width: 800px; height: %ipx">&nbsp;</div>', $chart->{id}."_".$chart->{type}, $height);

    print $fh
      sprintf('<div id="%s_table" style="width: 800px; height: 200px">&nbsp;</div>', $chart->{id})
      unless $chart->{no_table};
  }
  
  print $fh '</div>';
  print $fh "\n</div></body>\n</html>\n";

  return 1;
}


=head2 stats_html_head

  Arg 1      : arrayref $chart_data
  Example    : $html = $stats->stats_html_head($chart_data);
  Description: Generates the HTML header for HTML stats dump, includes
               JavaScript necessary for rendering and interacting with
               charts.
  Returntype : string
  Exceptions : none
  Caller     : dump_html()
  Status     : Stable

=cut

sub stats_html_head {
  my $self = shift;
  my $charts = shift;
  
  my ($js);
  foreach my $chart(@$charts) {
    my @keys = @{$self->sort_keys($chart->{data}, $chart->{sort})};
    
    my $type = ucfirst($chart->{type});
    
    # add colour
    if(defined($chart->{colours})) {
      my $co = 'slices: ['.join(", ", map { $chart->{colours}->{$_} ? '{color: "'.$chart->{colours}->{$_}.'"}' : '{}' } @keys).']';
      
      if(defined($chart->{options})) {
        $chart->{options} =~ s/}$/, $co}/;
      }
      else {
        $chart->{options} = "{$co}";
      }
    }
    
    # code to draw chart
    $js .= sprintf(
      "var %s = draw$type('%s', '%s', google.visualization.arrayToDataTable([['%s','%s'],%s]), %s);\n",
      $chart->{id}.'_'.$chart->{type},
      $chart->{id}.'_'.$chart->{type},
      $chart->{title},
      $chart->{header}->[0], $chart->{header}->[1],
      join(",", map {"['".$_."',".$chart->{data}->{$_}."]"} @keys),
      $chart->{options} || 'null',
    );
    
    unless($chart->{no_table}) {
      
      # code to draw table
      $js .= sprintf(
        "var %s = drawTable('%s', '%s', google.visualization.arrayToDataTable([['%s','%s'],%s]));\n",
        $chart->{id}.'_table',
        $chart->{id}.'_table',
        $chart->{title},
        $chart->{header}->[0], $chart->{header}->[1],
        join(",", map {"['".$_."',".$chart->{data}->{$_}."]"} @keys)
      );
      
      # interaction between table/chart
      $js .= sprintf(
        qq{
          google.visualization.events.addListener(%s, 'select', function() {
            %s.setSelection(%s.getSelection());
          });
          google.visualization.events.addListener(%s, 'select', function() {
            %s.setSelection(%s.getSelection());
          });
        },
        $chart->{id}.'_'.$chart->{type},
        $chart->{id}.'_table',
        $chart->{id}.'_'.$chart->{type},
        $chart->{id}.'_table',
        $chart->{id}.'_'.$chart->{type},
        $chart->{id}.'_table',
      );
    }
  }
  
  my $html =<<SHTML;
<html>
<head>
  <title>VEP summary</title>
  <script type="text/javascript" src="http://www.google.com/jsapi"></script>
  <script type="text/javascript">
    google.load('visualization', '1', {packages: ['corechart','table']});
  </script>
  <script type="text/javascript">
    
    function init() {
      // charts
      $js
    }
    
    function drawPie(id, title, data, options) {    
      var pie = new google.visualization.PieChart(document.getElementById(id));
      pie.draw(data, options);
      return pie;
    }
    function drawBar(id, title, data, options) {
      var bar = new google.visualization.ColumnChart(document.getElementById(id));
      bar.draw(data, options);
      return bar;
    }
    function drawTable(id, title, data) {
      var table = new google.visualization.Table(document.getElementById(id));
      table.draw(data, null);
      return table;
    }
    function drawLine(id, title, data, options) {
      var line = new google.visualization.LineChart(document.getElementById(id));
      line.draw(data, options);
      return line;
    }
    function drawArea(id, title, data, options) {
      var area = new google.visualization.AreaChart(document.getElementById(id));
      area.draw(data, options);
      return area;
    }
    google.setOnLoadCallback(init);
  </script>
  
  
  <style type="text/css">
    body {
      font-family: arial, sans-serif;
      margin: 0px;
      padding: 0px;
    }
    
    a {color: #36b;}
    a.visited {color: #006;}
    
    .stats_table {
      margin: 5px;
    }
    
    td {
      padding: 5px;
    }

    th.gradient {
      height: auto;
    }

    .stats_table tr:nth-child(odd) {
      background-color: #f0f0f0;
    }
    
    h3 {
      color: #666;
    }
    
    .masthead {
      background-color: white;
      color: rgb(204, 221, 255);
      height: 80px;
      width: 100%;
      padding: 0px;
    }
    
    .main {
      padding: 10px;
    }
    
    .gradient {
      background: #333366; /* Old browsers */
      background: -moz-linear-gradient(left,  #333366 0%, #ffffff 100%); /* FF3.6+ */
      background: -webkit-gradient(linear, left top, right top, color-stop(0%,#333366), color-stop(100%,#ffffff)); /* Chrome,Safari4+ */
      background: -webkit-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Chrome10+,Safari5.1+ */
      background: -o-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Opera 11.10+ */
      background: -ms-linear-gradient(left,  #333366 0%,#ffffff 100%); /* IE10+ */
      background: linear-gradient(to right,  #333366 0%,#ffffff 100%); /* W3C */
      filter: progid:DXImageTransform.Microsoft.gradient( startColorstr='#333366', endColorstr='#ffffff',GradientType=1 ); /* IE6-9 */
      
      padding: 0px;
      height: 80px;
      width: 700px;
    }
    
    .main_content {
      margin-left: 300px;
    }
    
    .sidemenu {
      width: 260px;
      position: fixed;
      border-style: solid;
      border-width: 2px;
      border-color: rgb(51, 51, 102);
    }
    
    .sidemenu_head {
      width: 250px;
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      padding: 5px;
    }
    
    .sidemenu_body {
      width: 250px;
      padding: 5px;
    }
  </style>
</head>
<body>
<div id="masthead" class="masthead">
  <div style="float: left; display: inline; padding: 10px; height: 80px;">
    <a href="http://www.ensembl.org/"><img src="http://static.ensembl.org/i/e-ensembl.png" width="185" height="45" alt="Ensembl"></a>
  </div>
  
  <div style="float: right; display: inline; height: 80px; background: white; padding: 10px;">
    <a href="http://www.ensembl.org/vep"><img src="http://www.ensembl.org/img/vep_logo.png" width="111" height="50" alt="VEP"></a>
  </div>
  <div class="gradient">
  </div>
</div>
<div class="main">
SHTML

  return $html;
}

1;
