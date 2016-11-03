#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

convert_cache.pl - a script to convert VEP caches to use tabix

http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#convert

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Storable qw(nstore_fd fd_retrieve freeze thaw);
use MIME::Base64;

use FindBin qw($RealBin);
use lib $RealBin;
use lib $RealBin.'/modules';

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = configure(scalar @ARGV);

# run the main sub routine
main($config);

sub configure {
  my $args = shift;
  
  my $config = {};
  
  GetOptions(
    $config,
    'help|h',            # displays help message
    'quiet|q',           # no output on STDOUT
    'force_overwrite|f', # force overwrite existing files
    'remove|r',          # remove cache files after we finish
    
    'species|s=s',       # species
    'dir|d=s',           # cache dir
    'version|v=s',       # version number
    
    'compress|c=s',      # eg zcat
    'bgzip|b=s',         # path to bgzip
    'tabix|t=s',         # path to tabix

    'sereal',            # convert transcript and reg caches to sereal
  ) or die "ERROR: Failed to parse command-line flags\n";
  
  # print usage message if requested or no args supplied
  if(defined($config->{help}) || !$args) {
    &usage;
    exit(0);
  }

  # check sereal
  if($config->{sereal}) {
    eval q{ use Sereal::Encoder; };
    die("ERROR: Could not use Sereal::Encoder perl module; perhaps you forgot to install it?\n$@") if $@;
  }
  
  $config->{dir} ||= join '/', ($ENV{'HOME'}, '.vep');
  die("ERROR: directory ".$config->{dir}." not found\n") unless -d $config->{dir};
  
  # check species
  opendir DIR, $config->{dir};
  my @species = grep {-d $config->{dir}.'/'.$_ && !/^\./} readdir DIR;
  closedir DIR;
  
  if(!defined($config->{species})) {
    my $msg = @species ? " or select one of the following:\n".join("\n", @species)."\n" : "";
    die("ERROR: No species specified (--species). Use \"--species all\"$msg\n");
  }
  elsif($config->{species} eq 'all') {
    $config->{species} = \@species;
  }
  else {
    $config->{species} = [split /\,/, $config->{species}];
    
    # check they exist
    foreach my $sp(@{$config->{species}}) {
      die("ERROR: Species $sp not found\n") unless grep {$sp eq $_} @species;
    }
  }
  
  # check versions
  my %versions;
  foreach my $sp(@{$config->{species}}) {
    opendir DIR, $config->{dir}.'/'.$sp;
    %{$versions{$sp}} = map {$_ => 1} grep {-d $config->{dir}.'/'.$sp.'/'.$_ && !/^\./} readdir DIR;
    closedir DIR;
  }
  
  if(!defined($config->{version})) {
    my $msg = keys %versions ? " or select one of the following:\n".join("\n", keys %versions)."\n" : "";
    die("ERROR: No version specified (--version). Use \"--version all\"$msg\n");
  }
  elsif($config->{version} eq 'all') {
    $config->{version} = \%versions;
  }
  else {
    $config->{version} = [split(/\,/, $config->{version})];
    
    # check they exist
    foreach my $v(@{$config->{version}}) {
      die("ERROR: Version $v not found\n") unless grep {defined($versions{$_}->{$v})} @{$config->{species}};
    }
    
    my %tmp;
    for my $sp(@{$config->{species}}) {
      %{$tmp{$sp}} = map {$_ => 1} @{$config->{version}};
    }
    $config->{version} = \%tmp;
  }
  
  $config->{compress} ||= 'gzip -dc';

  foreach my $tool(qw(bgzip tabix)) {
    unless($config->{$tool}) {
      if(`which $tool` =~ /$tool/) {
        $config->{$tool} = $tool;
      }
      elsif(-e $RealBin.'/htslib/'.$tool) {
        $config->{$tool} = $RealBin.'/htslib/'.$tool;
      }
      else {
        die("ERROR: Unable to convert cache without $tool\n");
      }
    }
  }
  
  return $config;
}

sub main {
  my $config = shift;
  
  my $bgzip    = $config->{bgzip};
  my $tabix    = $config->{tabix};
  my $zcat     = $config->{compress};
  my $base_dir = $config->{dir};
  
  my @files_to_remove;
  
  foreach my $sp(@{$config->{species}}) {
    debug($config, "Processing $sp");
    
    foreach my $v(keys %{$config->{version}->{$sp}}) {
      my $dir = join('/', ($base_dir, $sp, $v));
      $config->{dir} = $dir;
      
      next unless -d $dir;
      debug($config, "Processing version $v");
      
      opendir DIR, $dir;
      my @chrs = grep {-d $dir.'/'.$_ && !/^\./} readdir DIR;
      closedir DIR;
      
      # read cache info
      read_cache_info($config);
      
      # get pos col
      my @cols = @{$config->{cache_variation_cols}};
      my %var_cols = map {$cols[$_] => $_} (0..$#cols);
      $config->{pos_col} = $var_cols{start};
      
      foreach my $t($config->{sereal} ? qw(_tr _reg _var) : qw(_var)) {
      #foreach my $t(qw(_tr _reg _var)) {
        
        my %chr_files;
        my $total = 0;
        my $i = 0;
        
        debug($config, "Processing $t cache type");
        
        my $type = $t;
        my $orig_type = $type;
        $type = '' if $type eq '_tr';
        
        my $method = 'process'.$orig_type;
        my $method_ref = \&$method;
        
        foreach my $chr(@chrs) {    
          opendir DIR, $dir.'/'.$chr;
          my @files = grep {-f $dir.'/'.$chr.'/'.$_ && /\d+$type\.gz$/} readdir DIR;
          closedir DIR;
          
          # make sure we process in chromosomal order
          my @tmp;
          
          foreach my $file(@files) {
            if($file =~ m/(\d+)\-\d+$type\.gz/) {
              push @tmp, {
                s => $1,
                f => $file,
              };
            }
            else {
              die("ERROR: Filename $file doesn't look right\n");
            }
          }
          
          @files = map {$_->{f}} sort {$a->{s} <=> $b->{s}} @tmp;
          
          $chr_files{$chr} = \@files;
          $total += scalar @files;
        }
        
        foreach my $chr(@chrs) {
          
          my $outfilepath = join('/', ($dir, $chr, "all".$orig_type."s"));
          
          # check if files exist
          foreach my $file($outfilepath.'.gz', $outfilepath.'.gz.tbi') {
            if(-e $file) {
              if(defined($config->{force_overwrite})) {
                unlink($file) or die("ERROR: Failed to delete file $file\n");
              }
              else {
                die("ERROR: File $file already exists - use --force_overwrite to overwrite\n");
              }
            }
          }
          
          my $out_fh = new FileHandle;

          if($type eq '_var') {
            $out_fh->open(">$outfilepath") or die("ERROR: Could not write to file $outfilepath\n");
          }
          
          foreach my $file(@{$chr_files{$chr}}) {
            progress($config, $i++, $total);
            
            my $infilepath = join('/', ($dir, $chr, $file));
            
            &$method_ref($config, $chr, $infilepath, $out_fh);
            
            push @files_to_remove, $infilepath if defined($config->{remove});
          }
          
          $out_fh->close();
          
          # sort
          if($type eq '_var') {

            # bgzip
            my $bgzipout = `$bgzip $outfilepath 2>&1`;
            die("ERROR: bgzip failed\n$bgzipout") if $bgzipout;
            
            # tabix
            my ($b, $e) = ($config->{pos_col} + 2, $config->{pos_col} + 2);
            my $tabixout = `$tabix -s 1 -b $b -e $e $outfilepath\.gz 2>&1`;
            die("ERROR: tabix failed\n$tabixout") if $tabixout;
          }
        }
        
        end_progress($config);
      }
      
      $config->{cache_var_type} = 'tabix';
      $config->{cache_serialiser_type} = 'sereal' if $config->{sereal};
      $config->{cache_variation_cols} = 'chr,'.join(",", @{$config->{cache_variation_cols}}) unless $config->{cache_variation_cols}->[0] eq 'chr';
      
      open OUT, ">".$config->{dir}.'/info.txt' or die("ERROR: Could not write to info.txt\n");
      print OUT "# CACHE UPDATED ".get_time()."\n";
      foreach my $param(grep {/^cache\_/} keys %$config) {
        my $val = $config->{$param};
        if(ref($val) eq 'ARRAY') {
          $val = join ",", @$val;
        }
        $param =~ s/^cache\_//;
        print OUT "$param\t$val\n";
      }
      
      my $version_data = get_version_data($config);
      print OUT "source\_$_\t".$version_data->{$_}."\n" for keys %$version_data;
      
      close OUT;
    }
  }
  
  if(scalar @files_to_remove) {
    debug($config, "Removing old cache files");
    unlink($_) for @files_to_remove;
  }
  
  debug($config, "All done!");
}

sub process_tr {
  my ($config, $chr, $infilepath, $out_fh) = @_;
  my $zcat = $config->{compress};
  
  open my $fh, $zcat." ".$infilepath." |" or die("ERROR: Could not read from file $infilepath\n");
  my $tr_cache;
  $tr_cache = fd_retrieve($fh);
  close $fh;
  
  $config->{encoder} ||= Sereal::Encoder->new({compress => 1});

  $infilepath =~ s/\.gz/\.sereal/;

  die("ERROR: File $infilepath already exists - use --force_overwrite to overwrite\n") if !defined($config->{force_overwrite}) && -e $infilepath;

  open OUT, ">".$infilepath or die("ERROR: Could not write to dump file $infilepath");
  print OUT $config->{encoder}->encode($tr_cache);
  close OUT;
}

sub process_reg {
  my ($config, $chr, $infilepath, $out_fh) = @_;
  my $zcat = $config->{compress};
  
  open my $fh, $zcat." ".$infilepath." |" or die("ERROR: Could not read from file $infilepath\n");
  my $rf_cache;
  $rf_cache = fd_retrieve($fh);
  close $fh;
  
  $config->{encoder} ||= Sereal::Encoder->new({compress => 1});

  $infilepath =~ s/\.gz/\.sereal/;

  die("ERROR: File $infilepath already exists - use --force_overwrite to overwrite\n") if !defined($config->{force_overwrite}) && -e $infilepath;

  open OUT, ">".$infilepath or die("ERROR: Could not write to dump file $infilepath");
  print OUT $config->{encoder}->encode($rf_cache);
  close OUT;
}

sub process_var {
  my ($config, $chr, $infilepath, $out_fh) = @_;
  my $zcat = $config->{compress};
  
  my @tmp;
  
  open IN, "$zcat $infilepath | " or die("ERROR: Could not read from file $infilepath\n");
  while(<IN>) {
    chomp;
    my $l = length($_);
    
    while(1) {
      s/  / \. /g;
      s/ $/ \./g;
      last if length($_) == $l;
      $l = length($_);
    }
    
    s/ /\t/g;
    
    my @data = split /\t/, $_;
    my $pos = $data[$config->{pos_col}];
    
    push @tmp, {
      p => $pos,
      d => "$chr\t$_\n"
    }
  }
  close IN;
  
  print $out_fh $_ for map {$_->{d}} sort {$a->{p} <=> $b->{p}} @tmp;
}

# gets time
sub get_time() {
  my @time = localtime(time());

  # increment the month (Jan = 0)
  $time[4]++;

  # add leading zeroes as required
  for my $i(0..4) {
    $time[$i] = "0".$time[$i] if $time[$i] < 10;
  }

  # put the components together in a string
  my $time =
    ($time[5] + 1900)."-".
    $time[4]."-".
    $time[3]." ".
    $time[2].":".
    $time[1].":".
    $time[0];

  return $time;
}

# prints debug output with time
sub debug {
  my $config = shift;
  return if defined($config->{quiet});
  
  my $text = (@_ ? (join "", @_) : "No message");
  my $time = get_time;
  
  print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

# update or initiate progress bar
sub progress {
  my ($config, $i, $total) = @_;
  
  return if defined($config->{quiet}) || defined($config->{no_progress});
  
  $i = $total if $i > $total;
  
  my $width = $config->{terminal_width} || 60;
  my $percent = int(($i/$total) * 100);
  my $numblobs = int((($i/$total) * $width) - 2);
  
  # this ensures we're not writing to the terminal too much
  return if(defined($config->{prev_prog})) && $numblobs.'-'.$percent eq $config->{prev_prog};
  $config->{prev_prog} = $numblobs.'-'.$percent;
  
  #printf("\r%s of %s", $i, $total);
  printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

# end progress bar
sub end_progress {
  my $config = shift;
  return if defined($config->{quiet}) || defined($config->{no_progress});
  progress($config, 1,1);
  print "\n";
  delete $config->{prev_prog};
}

sub usage {
  print qq{#---------------#
# convert_cache.pl #
#---------------#

By Will McLaren (wm2\@ebi.ac.uk)

http://www.ensembl.org/info/docs/variation/vep/vep_script.html#convert_cache

Usage:
perl convert_cache.pl [arguments]
  
--help               -h   Print usage message and exit
--quiet              -q   Shhh!
--force_overwrite    -f   Overwrite existing cache files if found
--remove             -r   Remove old cache files after conversion

--dir [dir]          -d   Cache directory (default: \$HOME/.vep)
--species [species]  -s   Species cache to convert ("all" to do all found)
--version [version]  -v   Cache version to convert ("all" to do all found)

--compress [cmd]     -c   Path to binary/command to decompress gzipped files.
                          Defaults to "zcat", some systems may prefer "gzip -dc"
--bgzip [cmd]        -b   Path to bgzip binary (default: bgzip)
--tabix [cmd]        -t   Path to tabix binary (default: tabix)
};
}
