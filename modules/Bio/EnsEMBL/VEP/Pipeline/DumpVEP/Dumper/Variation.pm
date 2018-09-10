=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Variation;

use strict;
use warnings;

use FileHandle;

use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(trim_sequences get_matched_variant_alleles);
use Scalar::Util qw(looks_like_number);

our $CAN_USE_TABIX_PM;

BEGIN {
  if (eval { require Bio::EnsEMBL::IO::Parser::VCF4Tabix; 1 }) {
    $CAN_USE_TABIX_PM = 1;
  }
  else {
    $CAN_USE_TABIX_PM = 0;
  }
}

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper);

sub run {
  my $self = shift;

  my $vep_params = $self->get_vep_params();

  # make sure to include failed variants!
  $vep_params->{failed} = 1;

  my $config = Bio::EnsEMBL::VEP::Config->new($vep_params);

  my $region_size = $self->param('region_size');

  my $hive_dbc = $self->dbc;
  $hive_dbc->disconnect_if_idle() if defined $hive_dbc;

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
    config => $config,
    cache_region_size => $region_size,
  });

  my $cache = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation->new({
    config => $config,
    cache_region_size => $region_size,
    dir => $self->get_cache_dir($vep_params)
  });

  # create prefixed names here to use later
  if($vep_params->{freq_vcf} && !$self->{freq_vcf}) {
    foreach my $vcf_conf(@{$vep_params->{freq_vcf}}) {
      if(my $prefix = $vcf_conf->{prefix}) {
        $vcf_conf->{prefixed_pops} = [
          map {$_ ? $prefix.'_'.$_ : $prefix}
          @{$vcf_conf->{pops}}
        ];
      }
    }
    $self->{freq_vcf} = $vep_params->{freq_vcf};
  }

  # precache pubmed data
  $self->pubmed($as);

  $self->dump_chrs($as, $cache);

  # bgzip and tabix-index all_vars files
  if($self->param('convert')) {
    my $root = $self->get_cache_dir($vep_params);

    # need to find which column is start
    my @cols = @{$as->get_cache_columns()};
    my %indexes = map {$cols[$_] => $_} 0..$#cols;
    die("Cannot determine start column for indexing\n") unless defined($indexes{start});

    # add 2, 1 to correct for 0->1 indexing, 1 because we've added a chr column
    my $start_i = $indexes{start} + 2;

    foreach my $chr(keys %{{map {$_->{chr} => 1} @{$self->param('regions')}}}) {

      next unless -e "$root/$chr/all_vars";

      $self->run_system_command(
        sprintf(
          "cat %s/%s/all_vars | sort -k%i,%in | bgzip -c > %s/%s/all_vars.gz",
          $root, $chr, $start_i, $start_i, $root, $chr
        )
      );

      $self->run_system_command("tabix -C -s 1 -b $start_i -e $start_i $root/$chr/all_vars.gz");

      unlink("$root/$chr/all_vars");
    }
  }

  $self->dump_info($as, $self->get_cache_dir($vep_params));
  
  return;
}

sub dump_info {
  my ($self, $as, $dir) = @_;

  my $info_file = $dir.'/info.txt_variation';
  return if -e $info_file;

  my $fh = FileHandle->new;
  $fh->open(">$info_file");

  $self->_generic_dump_info($fh, $as, 0);

  $fh->close();

  $self->dump_info_converted($as, $dir) if $self->param('convert');
}

sub dump_info_converted {
  my ($self, $as, $dir) = @_;

  my $info_file = $dir.'/info.txt_variation_converted';
  return if -e $info_file;

  my $fh = FileHandle->new;
  $fh->open(">$info_file");

  $self->_generic_dump_info($fh, $as, 1);

  print $fh "var_type\ttabix\n";

  $fh->close();
}

sub _generic_dump_info {
  my ($self, $fh, $as, $converted) = @_;

  # var cache cols
  my @cols = (
    @{$as->get_cache_columns()},
    'pubmed',
    map {@{$_->{prefixed_pops} || $_->{pops}}} @{$self->{freq_vcf} || []}
  );
  unshift @cols, 'chr' if $converted;
  
  print $fh "variation_cols\t".join(",", @cols)."\n";

  my $info = $as->info;
  print $fh "source_$_\t".$info->{$_}."\n" for keys %$info;

  printf $fh "source_%s\t%s\n", $_->{name}, $_->{version} for
    grep {defined($_->{name}) && defined($_->{version})}
    @{$self->{freq_vcf} || []};
}

sub get_dumpable_object {
  my ($self, $as, $sr, $chr, $s) = @_;
  return $as->get_features_by_regions_uncached([[$sr, $s]], 1);
}

sub dump_obj {
  my ($self, $obj, $file, $chr) = @_;

  open DUMP, "| gzip -9 -c > ".$file or die "ERROR: Could not write to dump file $file\n";

  my $all_vars_fh;
  if($self->param('convert')) {
    my $all_vars_file = $file;
    $all_vars_file =~ s/[^\/]+?$//;
    $all_vars_file .= 'all_vars';

    $all_vars_fh = FileHandle->new();
    $all_vars_fh->open(">>$all_vars_file") or die $!;
  }

  # get freqs from VCFs?
  $self->freqs_from_vcf($obj, $chr) if $self->{freq_vcf};

  my $pubmed = $self->pubmed;

  foreach my $v(@$obj) {
    my @tmp = (
      $v->{variation_name},
      $v->{failed} == 0 ? '' : $v->{failed},
      $v->{somatic} == 0 ? '' : $v->{somatic},
      $v->{start},
      $v->{end} == $v->{start} ? '' : $v->{end},
      $v->{allele_string},
      $v->{strand} == 1 ? '' : $v->{strand},
      $v->{minor_allele} || '',
      defined($v->{minor_allele_freq}) && $v->{minor_allele_freq} =~ /^[0-9\.]+$/ ? sprintf("%.4f", $v->{minor_allele_freq}) : '',
      $v->{clin_sig} || '',
      $v->{phenotype_or_disease} == 0 ? '' : $v->{phenotype_or_disease},
    );

    push @tmp, $pubmed->{$v->{variation_name}} || '';

    if($self->{freq_vcf}) {
      foreach my $pop(map {@{$_->{prefixed_pops} || $_->{pops}}} @{$self->{freq_vcf}}) {
        push @tmp, $v->{$pop} || '';
      }
    }

    print DUMP join(" ", @tmp);
    print DUMP "\n";

    if($all_vars_fh) {
      print $all_vars_fh join("\t", $chr, map {$_ eq '' ? '.' : $_} @tmp);
      print $all_vars_fh "\n";
    }
  }

  close DUMP;
}

sub freqs_from_vcf {
  my $self = shift;
  my $obj = shift;
  my $chr = shift;

  die("ERROR: Cannot use freqs_from_vcf without Bio::EnsEMBL::IO::Parser::VCF4Tabix module\n") unless $CAN_USE_TABIX_PM;

  # sort by pos and exclude somatic vars
  my @list =
    grep {!$_->{somatic}}
    sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}}
    @$obj;
  return unless scalar @list;

  # put into a hash so we can lookup by pos
  my %by_pos;
  for(@list) {
    push @{$by_pos{$_->{start}}}, $_;
  }
  
  # iterate over each VCF file in the config
  foreach my $vcf_conf(@{$self->{freq_vcf}}) {
    my $file = $vcf_conf->{file};
    next unless -e $file;

    my $prefix = $vcf_conf->{prefix} || '';
    $prefix .= '_' if $prefix && $prefix !~ /\_$/;

    my $parser = $self->{_vcf_parsers}->{$file} ||= Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file);
    next unless $parser;

    $parser->seek($chr, $list[0]->{start} - 1, $list[-1]->{end} + 1);

    while($parser->next) {

      my $vcf_ref  = $parser->get_reference;
      my $vcf_pos  = $parser->get_raw_start;
      my @vcf_alts = @{$parser->get_alternatives};

      # scan from from pos to inferred end
      for my $start(grep {$by_pos{$_}} ($vcf_pos..($vcf_pos + length($vcf_ref)))) {

        foreach my $v(@{$by_pos{$start}}) {
          $DB::single = 1 if $v->{variation_name} eq 'TMP_ESP_1_179086420_179086420';

          my $matches = [];
          eval{
            $matches = get_matched_variant_alleles(
              {
                allele_string => $v->{allele_string},
                pos           => $v->{start},
                strand        => $v->{strand},
              },
              {
                ref  => $vcf_ref,
                alts => \@vcf_alts,
                pos  => $vcf_pos,
              }
            );
          };
          die "Failed to get vf matches for " .$v->{variation_name} ."\n" unless $@ eq '';  

          if(@$matches) {

            my %allele_map = map {$_->{b_index} => $_->{a_allele}} @$matches;

            # get INFO field data from VCF
            my $info = $parser->get_info;

            foreach my $pop(@{$vcf_conf->{pops}}) {

              my $info_prefix = '';
              my $info_suffix = '';

              # have to process ExAC differently from 1KG and ESP
              if($prefix =~ /exac|gnomad/i && $pop) {
                $info_suffix = '_'.$pop if $pop;
              }
              elsif($pop) {
                $info_prefix = $pop.'_';
              }

              my $tmp_f;

              if(exists($info->{$info_prefix.'AF'.$info_suffix})) {
                my $f = $info->{$info_prefix.'AF'.$info_suffix};
                my @split = split(',', $f);

                # there will be one item in @split for each of the original alts
                # since we may not be dealing with all the alts here
                # we have to use the indexes and alts we logged in %allele_map
                $tmp_f = join(',',
                  map {$allele_map{$_}.':'.($split[$_] == 0 ? 0 : sprintf('%.4g', $split[$_]))}
                  grep {$allele_map{$_}}
                  grep {looks_like_number($split[$_])}
                  0..$#split
                );
              }
              elsif(exists($info->{$info_prefix.'AC'.$info_suffix})) {
                my $c = $info->{$info_prefix.'AC'.$info_suffix};
                my $n = $info->{$info_prefix.'AN'.$info_suffix};
                my @split = split(',', $c);

                unless($n) {
                  $n += $_ for @split;
                }
                
                next unless $n;

                # ESP VCFs include REF as last allele, just to annoy everyone
                pop @split if scalar @split > scalar @vcf_alts;

                $tmp_f = join(',',
                  map {$allele_map{$_}.':'.($split[$_] ? sprintf('%.4g', $split[$_] / $n) : 0)}
                  grep {$allele_map{$_}}
                  grep {looks_like_number($split[$_])}
                  0..$#split
                );
              }              

              if(defined($tmp_f) && $tmp_f ne '') {
                my $store_name = $prefix.$pop;
                $store_name =~ s/\_$//;
                
                $v->{$store_name} = $v->{$store_name} ? $v->{$store_name}.','.$tmp_f : $tmp_f;
              }
            }
          }
        }
      }
    }
  }
}

# we want pubmed IDs by rsID
# this routine gets them from the DB in the first instance
# then caches this fetch to disk so future processes can read it
# as the query takes a little while to run
sub pubmed {
  my $self = shift;
  my $as = shift;

  if(!exists($self->{_pubmed})) {

    my %pm;
    my $file = sprintf(
      '%s/pubmed_%s_%s_%s.txt',
      $self->required_param('pipeline_dump_dir'),
      $self->required_param('species'),
      $self->required_param('ensembl_release'),
      $self->required_param('assembly')
    );
    my $lock = $file.'.lock';

    my $sleep_count = 0;
    if(-e $lock) {
      while(-e $lock) {
        sleep 1;
        die("I've been waiting for $lock to be removed for $sleep_count seconds, something may have gone wrong\n") if ++$sleep_count > 900;
      }
    }

    if(-e $file) {
      open IN, $file;
      while(<IN>) {
        chomp;
        my @split = split;
        $pm{$split[0]} = $split[1];
      }
      close IN;
    }

    elsif($as) {
      open OUT, ">$lock";
      print OUT "$$\n";
      close OUT;
      $self->{_lock_file} = $lock;

      my $sth = $as->get_adaptor('variation', 'variation')->dbc->prepare(qq{
        SELECT v.name, GROUP_CONCAT(p.pmid)
        FROM variation v, variation_citation c, publication p
        WHERE v.variation_id = c.variation_id
        AND c.publication_id = p.publication_id
        AND p.pmid IS NOT NULL
        GROUP BY v.variation_id
      });
      $sth->execute;

      my ($v, $p);
      $sth->bind_columns(\$v, \$p);

      open OUT, ">$file";

      while($sth->fetch()) {
        $pm{$v} = $p;
        print OUT "$v\t$p\n";
      }

      close OUT;

      unlink($lock);

      $sth->finish();
    }

    $self->{_pubmed} = \%pm;
  }

  return $self->{_pubmed};
}

sub DESTROY {
  unlink($_[0]->{_lock_file}) if $_[0]->{_lock_file};
}

1;
