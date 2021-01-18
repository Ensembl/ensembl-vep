#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

convert_cache.pl - a script to convert VEP caches to use tabix

http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#convert

=cut

use strict;
use Getopt::Long;
use FileHandle;
use File::Copy;
use IO::Socket;
use IO::Select;
use Storable qw(nstore_fd fd_retrieve freeze thaw);
use MIME::Base64;

use FindBin qw($RealBin);
use lib $RealBin;
use lib $RealBin.'/modules';

use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::CacheDir;

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = configure(scalar @ARGV);

# run the main sub routine
main($config);

sub configure {
  my $args = shift;
  
  my $config = {
    compress => 'gzip -dc',
    fork     => 1,
  };
  
  GetOptions(
    $config,
    'help|h',            # displays help message
    'quiet|q',           # no output on STDOUT
    'force_overwrite|f', # force overwrite existing files
    'remove|r',          # remove cache files after we finish
    'fork=i',            # fork
    
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
  my $version_count = 0;
  foreach my $sp(@{$config->{species}}) {
    opendir DIR, $config->{dir}.'/'.$sp;
    %{$versions{$sp}} = map {$_ => 1} grep {-d $config->{dir}.'/'.$sp.'/'.$_ && !/^\./} readdir DIR;
    closedir DIR;
    $version_count += keys %{$versions{$sp}};
  }

  die "ERROR: No valid directories found\n" unless $version_count;
  
  if(!defined($config->{version}) && $version_count > 1) {
    my $msg = keys %versions ? " or select one of the following:\n".join("\n",
      keys %{{
        map {$_ => 1}
        map {keys %{$versions{$_}}}
        keys %versions
      }}
    )."\n" : "";
    die("ERROR: No version specified (--version). Use \"--version all\"$msg\n");
  }
  elsif($config->{version} eq 'all' || $version_count == 1) {
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

  my $base_dir    = $config->{dir};
  my $fork_number = $config->{fork};
  
  foreach my $sp(@{$config->{species}}) {
    debug($config, "Processing $sp");
    
    foreach my $v(keys %{$config->{version}->{$sp}}) {
      my $dir = join('/', ($base_dir, $sp, $v));
      $config->{dir} = $dir;
      next unless -d $dir;

      debug($config, "Processing version $v");

      # read cache info
      my $config_obj = Bio::EnsEMBL::VEP::Config->new({dir => $dir, offline => 1, species => $sp});
      my $cache_dir_obj = Bio::EnsEMBL::VEP::CacheDir->new({dir => $dir, config => $config_obj});

      # work out which types we're processing based on cache info and user options
      my @types;
      push @types, '_var' unless ($cache_dir_obj->info->{var_type} || '') eq 'tabix' || !$cache_dir_obj->info->{variation_cols};
      if($config->{sereal}) {
        push @types, qw(_tr _reg) unless ($cache_dir_obj->info->{serialiser_type} || '') eq 'sereal';
      }

      unless(@types) {
        debug($config, "No unprocessed types remaining, skipping");
        next;
      }

      opendir DIR, $dir;
      my @chrs = grep {-d $dir.'/'.$_ && !/^\./} readdir DIR;
      closedir DIR;
      
      # get pos col
      my @cols = @{$cache_dir_obj->info->{variation_cols}};
      my %var_cols = map {$cols[$_] => $_} (0..$#cols);
      $config->{pos_col} = $var_cols{start};
      
      foreach my $t(@types) {
        my %chr_files;
        my $total = 0;
        my $i = 0;
        
        debug($config, "Processing $t cache type");

        my $chr_files = get_chr_files($dir, \@chrs, $t);
        my $total = 0;
        my $i = 0;
        $total += scalar @{$chr_files->{$_}} for keys %$chr_files;
        
        if($fork_number <= 1) {

          foreach my $chr(@chrs) {
            progress($config, $i, $total);
            process_chr_type($config, $dir, $chr, $t, $chr_files->{$chr});
            $i += scalar @{$chr_files->{$chr}};
          }
        }
        else {
          my $active_forks = 0;
          my @pids;
          my $sel = IO::Select->new;

          while(keys %$chr_files or $active_forks) {

            # only spawn new forks if we have space
            if($active_forks <= $fork_number) {

              my $chr = (keys %$chr_files)[0];
              my $files = delete($chr_files->{$chr});

              # create sockets for IPC
              my ($child, $parent);
              socketpair($child, $parent, AF_UNIX, SOCK_STREAM, PF_UNSPEC) or throw("ERROR: Failed to open socketpair: $!");
              $child->autoflush(1);
              $parent->autoflush(1);
              $sel->add($child);

              # fork
              my $pid = fork;
              if(!defined($pid)) {
                throw("ERROR: Failed to fork\n");
              }
              elsif($pid) {
                push @pids, $pid;
                $active_forks++;
              }
              elsif($pid == 0) {
                process_chr_type($config, $dir, $chr, $t, $files);
                print $parent scalar @$files;
                exit(0);
              }
            }

            # read child input
            while(my @ready = $sel->can_read()) {
              my $no_read = 1;

              foreach my $fh(@ready) {
                $no_read++;

                my $line = join('', $fh->getlines());
                next unless $line;
                $no_read = 0;

                # finish up
                $sel->remove($fh);
                $fh->close;
                $active_forks--;
                progress($config, $i, $total);
                $i += $line;
              }

              # read-through detected, DIE
              die("ERROR: Forked process(es) died\n") if $no_read;

              last if $active_forks < $fork_number;
            }
          }

          waitpid($_, 0) for @pids;
        }
        
        end_progress($config);
      }

      open OUT, ">".$config->{dir}.'/info.txt.new' or die("ERROR: Could not write to info.txt.new\n");
      print OUT "# CACHE UPDATED ".get_time()."\n";

      open IN, $config->{dir}.'/info.txt' or die $!;
      while(<IN>) {
        if(/^variation_cols/) {
          print OUT "variation_cols\t".join(',', 'chr', @{$cache_dir_obj->info->{variation_cols}})."\n";
        }
        elsif(/\# CACHE UPDATED/) {
          next;
        }
        else {
          print OUT $_;
        }
      }

      print OUT "var_type\ttabix\n";
      print OUT "serialiser_type\tsereal\n" if $config->{sereal};

      close OUT;
      close IN;

      move($config->{dir}.'/info.txt', $config->{dir}.'/info.txt.bak');
      move($config->{dir}.'/info.txt.new', $config->{dir}.'/info.txt');
    }
  }
  
  debug($config, "All done!");
}

sub get_chr_files {
  my ($dir, $chrs, $type) = @_;

  my %chr_files;
  $type = '' if $type eq '_tr';

  foreach my $chr(@$chrs) {
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
  }

  return \%chr_files;
}

sub process_chr_type {
  my ($config, $dir, $chr, $type, $files) = @_;
  
  my $bgzip = $config->{bgzip};
  my $tabix = $config->{tabix};
  my $zcat  = $config->{compress};

  my $orig_type = $type;
  $type = '' if $type eq '_tr';
  
  my $method = 'process'.$orig_type;
  my $method_ref = \&$method;
          
  my $outfilepath = join('/', ($dir, $chr, "all".$orig_type."s"));

  my @files_to_remove;
  
  # check if files exist
  foreach my $file($outfilepath.'.gz', $outfilepath.'.gz.tbi', $outfilepath.'.gz.csi') {
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
  
  foreach my $file(@$files) {    
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
    my $tabixout = `$tabix -C -s 1 -b $b -e $e $outfilepath\.gz 2>&1`;
    die("ERROR: tabix failed\n$tabixout") if $tabixout;
  }
  
  if(scalar @files_to_remove) {
    unlink($_) for @files_to_remove;
  }
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

http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#convert

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
                          Defaults to "gzip -dc", some systems may prefer "zcat"
--bgzip [cmd]        -b   Path to bgzip binary (default: bgzip)
--tabix [cmd]        -t   Path to tabix binary (default: tabix)
};
}
