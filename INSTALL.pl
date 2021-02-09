#!/usr/bin/env perl

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

=head1 NAME

INSTALL.pl - a script to install required code and data for VEP

=cut

use strict;
use FindBin qw($RealBin);
use lib $RealBin.'/modules';
use Getopt::Long;
use File::Path qw(mkpath rmtree);
use File::Copy;
use File::Basename;
use Net::FTP;
use Cwd;
use Scalar::Util qw(looks_like_number);
use Bio::EnsEMBL::VEP::Utils qw(get_version_data get_version_string);

our (
  $DEST_DIR,
  $ENS_CVS_ROOT,
  $API_VERSION,
  $DATA_VERSION,
  $ASSEMBLY,
  $ENS_GIT_ROOT,
  $BIOPERL_URL,
  $CACHE_URL,
  $CACHE_URL_INDEXED,
  $CACHE_DIR,
  $PLUGINS,
  $PLUGIN_URL,
  $PLUGINS_DIR,
  $FASTA_URL,
  $FTP_USER,
  $HELP,
  $NO_UPDATE,
  $SPECIES,
  $AUTO,
  $QUIET,
  $PREFER_BIN,
  $CONVERT,
  $TEST,
  $NO_HTSLIB,
  $LIB_DIR,
  $HTSLIB_DIR,
  $BIODBHTS_DIR,
  $REALPATH_DEST_DIR,
  $NO_TEST,
  $NO_BIOPERL,
  $ua,

  $CAN_USE_CURL,
  $CAN_USE_LWP,
  $CAN_USE_HTTP_TINY,
  $CAN_USE_ARCHIVE,
  $CAN_USE_UNZIP,
  $CAN_USE_GZIP,
  $CAN_USE_TAR,
  $CAN_USE_DBI,
  $CAN_USE_DBD_MYSQL,
);


## VERSIONS OF INSTALLED SOFTWARE
## MAY BE UPDATED IF SUCCESSFULLY TESTED
########################################
our $HTSLIB_VERSION  = '1.9';             # latest release as of release/98
our $BIOHTS_VERSION  = '2.11';            # latest 2.X release as of release/98
our $BIOPERL_VERSION = 'release-1-6-924'; # frozen, no pressing need to update


## BEGIN BLOCK, CHECK WHAT MODULES ETC WE CAN USE
#################################################

BEGIN {
  if(eval q{ use LWP::Simple qw(getstore get $ua); 1 }) {
    $CAN_USE_LWP = 1;

    # set up a user agent's proxy (excluding github)
    $ua->env_proxy;
  }

  if(eval q{ use DBI; 1 }) {
    $CAN_USE_DBI = 1;
  }

  if(eval q{ use DBD::mysql; 1 }) {
    $CAN_USE_DBD_MYSQL = 1;
  }

  $CAN_USE_CURL      = 1 if `which curl` =~ /\/curl/;
  $CAN_USE_HTTP_TINY = 1 if eval q{ use HTTP::Tiny; 1 };
  $CAN_USE_ARCHIVE   = 1 if eval q{ use Archive::Extract; 1 };
  $CAN_USE_UNZIP     = 1 if `which unzip` =~ /\/unzip/;
  $CAN_USE_GZIP      = 1 if `which gzip` =~ /\/gzip/;
  $CAN_USE_TAR       = 1 if `which tar` =~ /\/tar/;
}

$| = 1;

# CONFIGURE
###########

# other global data
my @API_MODULES = (
  { name => 'ensembl',           path => ' ',         test_pm => 'Bio::EnsEMBL::Registry' },
  { name => 'ensembl-variation', path => 'Variation', test_pm => 'Bio::EnsEMBL::Variation::Variation' },
  { name => 'ensembl-funcgen',   path => 'Funcgen',   test_pm => 'Bio::EnsEMBL::Funcgen::RegulatoryFeature' },
  { name => 'ensembl-io',        path => 'IO,Utils',  test_pm => 'Bio::EnsEMBL::IO::Parser' },
);
my $ensembl_url_tail = '/archive/';
my $archive_type = '.zip';
my $git_api_root = 'https://api.github.com/repos/Ensembl/';
my $VEP_MODULE_NAME = 'ensembl-vep';

our (@store_species, @indexes, @files, $ftp, $dirname);

GetOptions(
  'DESTDIR|d=s'        => \$DEST_DIR,
  'VERSION|v=i'        => \$API_VERSION, # Deprecated
  'CACHE_VERSION|e=i'  => \$DATA_VERSION,
  'ASSEMBLY|y=s'       => \$ASSEMBLY,
  'BIOPERL|b=s'        => \$BIOPERL_URL,
  'CACHEURL|u=s'       => \$CACHE_URL,
  'CACHEDIR|c=s'       => \$CACHE_DIR,
  'FASTAURL|f=s'       => \$FASTA_URL,
  'HELP|h'             => \$HELP,
  'NO_UPDATE|n'        => \$NO_UPDATE,
  'SPECIES|s=s'        => \$SPECIES,
  'PLUGINS|g=s'        => \$PLUGINS,
  'PLUGINSDIR|r=s'     => \$PLUGINS_DIR,
  'PLUGINURL=s'        => \$PLUGIN_URL,
  'AUTO|a=s'           => \$AUTO,
  'QUIET|q'            => \$QUIET,
  'PREFER_BIN|p'       => \$PREFER_BIN,
  'CONVERT|t'          => \$CONVERT,
  'TEST'               => \$TEST,
  'NO_HTSLIB|l'        => \$NO_HTSLIB,
  'NO_TEST'            => \$NO_TEST,
  'NO_BIOPERL'         => \$NO_BIOPERL
) or die("ERROR: Failed to parse arguments");

# load version data
our $CURRENT_VERSION_DATA = get_version_data($RealBin.'/.version');
our $VERSION = $CURRENT_VERSION_DATA->{$VEP_MODULE_NAME}->{release};
$VERSION =~ s/release\///;

if($HELP) {
  usage();
  exit(0);
}

# check user has DBI and DBD::mysql
die(
  "ERROR: DBI module not found. VEP requires the DBI perl module to function\n\n".
  "http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#requirements\n"
) unless $CAN_USE_DBI;

warn(
  "WARNING: DBD::mysql module not found. VEP can only run in offline (--offline) mode without DBD::mysql installed\n\n".
  "http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#requirements\n"
) unless $CAN_USE_DBD_MYSQL;


my $default_dir_used = check_default_dir();

$LIB_DIR            = $DEST_DIR;
$HTSLIB_DIR         = $LIB_DIR.'/htslib';
$BIODBHTS_DIR       = $LIB_DIR.'/biodbhts';
$REALPATH_DEST_DIR .= Cwd::realpath($DEST_DIR).'/Bio';
$DEST_DIR          .= '/Bio';
$dirname            = dirname(__FILE__) || '.';

$ENS_GIT_ROOT ||= 'https://github.com/Ensembl/';
$BIOPERL_URL  ||= "https://github.com/bioperl/bioperl-live/archive/$BIOPERL_VERSION.zip";
$API_VERSION  ||= $CURRENT_VERSION_DATA->{$VEP_MODULE_NAME}->{release};
$DATA_VERSION ||= $API_VERSION;
$CACHE_DIR    ||= $ENV{HOME} ? $ENV{HOME}.'/.vep' : 'cache';
$PLUGINS_DIR  ||= $CACHE_DIR.'/Plugins';
$FTP_USER     ||= 'anonymous';

## Set the indexed cache url if it's been overwritten by the user
$CACHE_URL_INDEXED = $CACHE_URL;

$CACHE_URL  ||= "ftp://ftp.ensembl.org/pub/release-$DATA_VERSION/variation/vep";
$CACHE_URL_INDEXED  ||= "ftp://ftp.ensembl.org/pub/release-$DATA_VERSION/variation/indexed_vep_cache";
$FASTA_URL  ||= "ftp://ftp.ensembl.org/pub/release-$DATA_VERSION/fasta/";
$PLUGIN_URL ||= 'https://raw.githubusercontent.com/Ensembl/VEP_plugins';

# using PREFER_BIN can save memory when extracting archives
$PREFER_BIN = 0 unless defined($PREFER_BIN);
$Archive::Extract::PREFER_BIN = $PREFER_BIN == 0 ? 0 : 1;

$QUIET = 0 unless $AUTO;

# updates to ensembl-vep available?
update() unless $NO_UPDATE;

# auto?
if($AUTO) {

  # check
  die("ERROR: Failed to parse AUTO string - must contain any of a (API), l (FAIDX/htslib), c (cache), f (FASTA), p (plugins)\n") unless $AUTO =~ /^[alcfp]+$/i;

  # require species
  if($AUTO =~ /[cf]/i) {
    die("ERROR: No species specified\n") unless $SPECIES;
    $SPECIES = [split /\,/, $SPECIES];
  }

  # require plugin list
  if($AUTO =~ /p/i) {
    die("ERROR: No plugins specified\n") unless $PLUGINS;
    $PLUGINS = [split /\,/, $PLUGINS];
  }

  # run subs
  if($AUTO =~ /l/ && $AUTO !~ /a/) {
    my $curdir = getcwd;
    chdir $curdir;
    install_biodbhts();
    chdir $curdir;

    # remove Bio dir if empty
    opendir DIR, $DEST_DIR;
    my @files = grep {!/^\./} readdir DIR;
    closedir DIR;

    if(scalar @files <= 1) {
      rmtree($DEST_DIR.'/'.$files[0]);
      rmtree($DEST_DIR);
    }
  }

  api()   if $AUTO =~ /a/;
  cache() if $AUTO =~ /c/;
  fasta() if $AUTO =~ /f/;
  plugins() if $AUTO =~ /p/;
}

else {
  print "\nHello! This installer is configured to install v$API_VERSION of the Ensembl API for use by the VEP.\nIt will not affect any existing installations of the Ensembl API that you may have.\n\nIt will also download and install cache files from Ensembl's FTP server.\n\n" unless $QUIET;

  # run subs
  api() if check_api();
  cache();
  fasta();
  plugins();
}


# clean up
if(-d "$CACHE_DIR/tmp" && !$TEST) {
  rmtree("$CACHE_DIR/tmp") or die "ERROR: Could not delete directory $CACHE_DIR/tmp\n";
}

print "\nAll done\n" unless $QUIET;


##########################################################################
##########################################################################
##########################################################################


# UPDATE
########
sub update() {

  my $module = $VEP_MODULE_NAME;

  # check for major version update
  my $repo_file = "$RealBin/$$.repo_file";
  download_to_file(
    "$git_api_root$module",
    $repo_file
  );

  my $default_branch;
  open IN, $repo_file;
  while(<IN>) {
    if(/default_branch.+\:.+\"(.+?)\"/) {
      $default_branch = $1;
      last;
    }
  }
  close IN;

  unlink($repo_file);

  unless($default_branch) {
    print "WARNING: Unable to carry out version check for '$module'\n" unless $QUIET;
    return;
  }

  my $default_branch_number = $default_branch;
  $default_branch_number =~ s/release\/// if $default_branch_number;

  my $current_branch = $CURRENT_VERSION_DATA->{'ensembl-vep'}->{release};

  # Check if the $API_VERSION has been set by the deprecated "--VERSION" flag
  my $api_branch  = $API_VERSION;
     $api_branch  =~ s/release\///;

  # branch provided by the "--VERSION" flag
  if ($api_branch != $current_branch) {
    print "The 'VERSION' installation flag has been deprecated.\n\n";
    my $branch = looks_like_number($API_VERSION) ? 'release/'.$API_VERSION : $API_VERSION;
    if(`which git` && -d $RealBin.'/.git') {
      print "Please, use git to update '$module' using the commands:\n\n";
      print "\tgit pull\n";
      print "\tgit checkout $branch\n\n";
    }
    else {
      print "Please, re-download '$module' with the desired version.\n\n";
    }
    print "Exit VEP install script\n";
    exit(0);
  }

  my $message;

  # don't have latest
  if($current_branch < $default_branch_number) {
    $message = 
      "Version check reports a newer release of '$module' is available ".
      "(installed: $current_branch, available: $default_branch_number)\n";
  }

  # do have latest, but there might be updates
  elsif($current_branch == $default_branch_number) {
    my $git_sub = get_vep_sub_version($current_branch);
    my $have_sub = $CURRENT_VERSION_DATA->{$module}->{sub};

    $message = sprintf(
      "Version check reports there are post-release updates available of %s (installed: %s.%.7s, available: %s.%.7s)\n",
      $module,
      $current_branch, $have_sub,
      $current_branch, $git_sub
    ) unless $git_sub eq $have_sub;
  }

  if($message) {
    print "\n$message\n";

    # user has git, suggest they use that instead
    if(`which git` && -d $RealBin.'/.git') {
      print "We recommend using git to update '$module', by exiting this installer and running:\n\n";
      print "\tgit pull\n";
      print "\tgit checkout $default_branch\n" if $current_branch ne $default_branch_number;
    }
    else {
      print "You should exit this installer and re-download '$module' if you wish to update\n";
    }

    print "\nDo you wish to exit so you can get updates (y) or continue (n): ";

    my $ok = <>;

    if($ok !~ /^n/i) {
      print "OK, bye!\n";
      print "\nNB: Remember to re-run INSTALL.pl after updating to check for API updates\n";
      exit(0);
    }
  }
  else {
    return;
  }
}


# CHECKS DIR SETUP AND PATHS ETC
################################
sub check_default_dir {
  my $this_os =  $^O;
  my $default_dir_used;
  # check if $DEST_DIR is default
  if(defined($DEST_DIR)) {
    print "Using non-default API installation directory $DEST_DIR.\n";
    print "Please note this just specifies the location for downloaded API files. The vep script will remain in its current location where ensembl-vep was unzipped.\n";
    if(!defined($AUTO)){
      print "Have you \n";
      print "1. added $DEST_DIR to your PERL5LIB environment variable?\n";
      print "2. added $DEST_DIR/htslib to your PATH environment variable?\n";
      if( $this_os eq 'darwin' && !$NO_HTSLIB) {
        print "3. added $DEST_DIR/htslib to your DYLD_LIBRARY_PATH environment variable?\n";
      }
      print "(y/n): ";

      my $ok = <>;
      if($ok !~ /^y/i) {
        print "Exiting. Please \n";
        print "1. add $DEST_DIR to your PERL5LIB environment variable\n";
        print "2. add $DEST_DIR/htslib to your PATH environment variable\n";
        if( $this_os eq 'darwin' && !$NO_HTSLIB) {
          print "3. add $DEST_DIR/htslib to your DYLD_LIBRARY_PATH environment variable\n";
        }
        exit(0);
      }
    }
    else {
      print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
      print "PLEASE REMEMBER TO \n";
      print "1. add $DEST_DIR to your PERL5LIB environment variable\n";
      print "2. add $DEST_DIR/htslib to your PATH environment variable\n";
      if( $this_os eq 'darwin' && !$NO_HTSLIB) {
        print "3. add $DEST_DIR/htslib to your DYLD_LIBRARY_PATH environment variable\n";
        print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
      }
    }
    if( ! -d $DEST_DIR ) {
      mkdir $DEST_DIR || die "Could not make destination directory $DEST_DIR"
    }
    $default_dir_used = 0;
  }

  else {
    $DEST_DIR ||= '.';
    $default_dir_used = 1;
    my $current_dir = cwd();

    if( !$NO_HTSLIB && $this_os eq 'darwin' ) {
      print "Installation on OSX requires that you set up some paths before running this installer.\n";
      if(!defined($AUTO)){
        print "Have you \n";
        print "1. added $current_dir/htslib to your DYLD_LIBRARY_PATH environment variable?\n";
        print "(y/n): ";
        my $ok = <>;
        if($ok !~ /^y/i) {
          print "Exiting. Please \n";
          print "1. add $current_dir/htslib to your DYLD_LIBRARY_PATH environment variable\n";
          exit(0);
        }
      }
      else{
        print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
        #print "\nPLEASE REMEMBER TO ADD $current_dir/htslib TO YOUR DYLD_LIBRARY_PATH ENVIRONMENT VARIABLE\n";
        print "\nPLEASE REMEMBER TO \n";
        print "1. add $current_dir/htslib to your DYLD_LIBRARY_PATH environment variable\n";
        print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
      }
    }
  }

  return $default_dir_used;
}


# API
#####
sub api() {
  setup_dirs();
  my $curdir = getcwd;
  unless($NO_BIOPERL) {
    bioperl();
  }  

  # htslib needs to find bioperl to pass tests
  $ENV{PERL5LIB} = $ENV{PERL5LIB} ? $ENV{PERL5LIB}.':'.$DEST_DIR : $DEST_DIR;

  unless($NO_HTSLIB) {
    chdir $curdir;
    install_biodbhts();
  }

  chdir $curdir;
  install_api();
  test() unless $NO_TEST;
}


# CHECK EXISTING
################
sub check_api() {
  print "Checking for installed versions of the Ensembl API..." unless $QUIET;

  my $has_api = {};
  my $updates = {};
  my $core_version;
  my @unknown_versions = ();

  foreach my $module_hash(@API_MODULES) {

    my $module  = $module_hash->{name};
    my $test_pm = $module_hash->{test_pm};

    eval "require $test_pm";
    $has_api->{$module} = $@ ? 0 : 1;
    
    if($has_api->{$module}) {
      my $have_sub = $CURRENT_VERSION_DATA->{$module} ? ($CURRENT_VERSION_DATA->{$module}->{sub} || '') : '';
      my $git_sub = get_module_sub_version($module);

      if($have_sub) {
        $updates->{$module} = [$have_sub, $git_sub] if $have_sub ne $git_sub;
      }
      else {
        push @unknown_versions, $module;
      }

      if($module eq 'ensembl') {
        $core_version = Bio::EnsEMBL::Registry->software_version;
      }
    }
  }

  print "done\n";

  my $total = 0;
  $total += $_ for values %$has_api;

  my $message;

  if($total == scalar @API_MODULES) {

    if(defined($core_version)) {
      if(!looks_like_number($API_VERSION)) {
        $message = "Your reported version ($API_VERSION) is a non-standard release number";
      }
      elsif($core_version == $API_VERSION) {

        if(scalar keys %$updates) {
          $message =
            "There are updates available for these modules:\n  ".
            join(
              "\n  ",
              map {
                sprintf(
                  "%-20s : installed = %.7s, available = %.7s",
                  $_, $updates->{$_}->[0], $updates->{$_}->[1]
                )
              } keys %$updates
            );
        }
        else {
          $message = "It looks like you already have v$API_VERSION of the API installed.\n";

          if(@unknown_versions) {
            $message .=
              "No version information was available for the following modules:\n".
              join("\n", map {" - ".$_} @unknown_versions).
              "\nUpdates may be available but you will need to perform these manually.";
          }
          else {
            $message .= "You shouldn't need to install the API";
          }
        }
      }

      elsif($core_version > $API_VERSION) {
        $message = "It looks like this installer is for an older distribution ($API_VERSION) of the API than you already have ($core_version)";
      }

      else {
        $message = "It looks like you have an older version ($core_version) of the API installed.\nThis installer will install a limited set of the API v$API_VERSION for use by the VEP only";
      }
    }

    else {
      $message = "It looks like you have an unidentified version of the API installed.\nThis installer will install a limited set of the API v$API_VERSION for use by the VEP only"
    }
  }

  elsif($total > 0) {
    $message = "It looks like you already have the following API modules installed:\n\n".(join "\n", grep {$has_api->{$_}} keys %$has_api)."\n\nThe VEP requires the ensembl, ensembl-io, ensembl-variation and ensembl-funcgen modules";
  }

  if($AUTO =~ /a/ || !defined($message)) {
    return 1;
  }
  else {
    print $message unless $QUIET;

    print "\n\nSkip to the next step (n) to install cache files\n\nDo you want to continue installing the API (y/n)? ";
    my $ok = <>;

    if($ok !~ /^y/i) {
      print " - skipping API installation\n" unless $QUIET;
      return 0;
    }
    else {
      return 1;
    }
  }
}


# SETUP
#######
sub setup_dirs() {

  print "\nSetting up directories\n" unless $QUIET;

  # check if install dir exists
  if(-e $DEST_DIR) {
    my $ok;

    if($AUTO) {
      $ok = 'y';
    }
    else {
      print "Destination directory $DEST_DIR already exists.\nDo you want to overwrite it (if updating VEP this is probably OK) (y/n)? ";

      $ok = <>;
    }

    if($ok !~ /^y/i) {
      print "Exiting\n";
      exit(0);
    }

    else {
      unless($default_dir_used || $AUTO) {
        print "WARNING: You are using a non-default install directory.\nPressing \"y\" again will remove $DEST_DIR and its contents!!!\nAre you really, really sure (y/n)? ";
        $ok = <>;

        if($ok !~ /^y/i) {
          print "Exiting\n";
          exit(0);
        }
      }

      # try to delete the existing dir
      rmtree($DEST_DIR) or die "ERROR: Could not delete directory $DEST_DIR\n";
    }
  }

  mkdir($DEST_DIR) or die "ERROR: Could not make directory $DEST_DIR\n";
  mkdir($DEST_DIR.'/tmp') or die "ERROR: Could not make directory $DEST_DIR/tmp\n";
}


# INSTALL API
#############
sub install_api() {

  print "\nDownloading required Ensembl API files\n" unless $QUIET;

  my $release_url_string = looks_like_number($API_VERSION) ? 'release/'.$API_VERSION : $API_VERSION;
  my $release_path_string = looks_like_number($API_VERSION) ? 'release-'.$API_VERSION : $API_VERSION;

  foreach my $module_hash(@API_MODULES) {
    my $module = $module_hash->{name};

    # do we need to update this?
    my $have_sub = $CURRENT_VERSION_DATA->{$module} ? ($CURRENT_VERSION_DATA->{$module}->{sub} || '') : '';

    my $url = $ENS_GIT_ROOT.$module.$ensembl_url_tail.$release_url_string.$archive_type;

    print " - fetching $module\n" unless $QUIET;
    my $target_file = $DEST_DIR.'/tmp/'.$module.$archive_type;
    mkdir($DEST_DIR.'/tmp/') unless -d $DEST_DIR.'/tmp/';
    download_to_file($url, $target_file) unless -e $target_file;

    print " - unpacking $target_file\n" unless $QUIET;
    unpack_arch("$DEST_DIR/tmp/$module$archive_type", "$DEST_DIR/tmp/");

    print " - moving files\n" unless $QUIET;
    foreach my $module_path (split(',',$module_hash->{path})) {
      my $module_dir_suffix = $module_path eq ' ' ? '' : '/'.$module_path;
      my $module_dir_from   = "$DEST_DIR/tmp/$module\-$release_path_string/modules/Bio/EnsEMBL$module_dir_suffix";
      my $module_dir_to     = "$DEST_DIR/EnsEMBL$module_dir_suffix";

      # If the target directory already exist, we can't overwrite it.
      if (-d $module_dir_to) {
        # One solution is to move the content of the directory instead.
        # However we need to loop over the files/directories within the module directory because the 'move()' method doesn't allow wildcards.
        opendir DH, $module_dir_from;
        while(my $file_or_dir = readdir DH) {
          next if ($file_or_dir =~ /^\.+$/);
          move("$module_dir_from/$file_or_dir", "$module_dir_to/$file_or_dir") or die "ERROR: Could not move '$module_dir_from/$file_or_dir'\n".$!;
        }
        closedir DH;
      }
      else {
        move($module_dir_from, $module_dir_to) or die "ERROR: Could not move the directory '$module_dir_from'\n".$!;
      }
    }

    # now get latest commit from github API
    print " - getting version information\n" unless $QUIET;
    my $git_sub = get_module_sub_version($module);

    mkdir("$RealBin/.version/") unless -d "$RealBin/.version/";
    open OUT, ">$RealBin/.version/$module" or die $!;
    print OUT "release $API_VERSION\nsub $git_sub\n";
    close OUT;

    rmtree("$DEST_DIR/tmp/$module\-$release_path_string") or die "ERROR: Failed to remove directory: $!\n";
  }
}

sub get_module_sub_version {
  my $module = shift;

  my $sub_file = "$RealBin/$$\.$module.sub";
  my $release_url_string = looks_like_number($API_VERSION) ? 'release/'.$API_VERSION : $API_VERSION;

  download_to_file(
    "$git_api_root$module/commits?sha=$release_url_string",
    $sub_file
  );

  open IN, $sub_file or die $!;
  my $sub;
  while(<IN>) {
    if(/\"sha\": \"(.+?)\"/) {
      $sub = $1;
      last;
    }
  }
  close IN;

  unlink($sub_file);

  return $sub;
}

sub get_vep_sub_version {
  my $release = shift || $API_VERSION;

  my $sub_file = "$RealBin/$$\.$VEP_MODULE_NAME.sub";
  my $release_url_string = looks_like_number($release) ? 'release/'.$release : $release;

  download_to_file(
    sprintf(
      'https://raw.githubusercontent.com/Ensembl/%s/%s/modules/Bio/EnsEMBL/VEP/Constants.pm',
      $VEP_MODULE_NAME,
      $release_url_string
    ),
    $sub_file
  );

  open IN, $sub_file or die $!;
  my $sub;
  while(<IN>) {
    if(/VEP_SUB_VERSION \= (.+)\;/) {
      $sub = $1;
      last;
    }
  }
  close IN;

  unlink($sub_file);

  return $sub;
}

# HTSLIB download/make
######################
sub install_htslib() {

  #actually decided to follow Bio::DB::Sam template
  # STEP 0: various dependencies
  my $git = `which git`;
  $git or die <<END;
  'git' command not in path. Please install git and try again.
  (or to skip Bio::DB::HTS/htslib install re-run with --NO_HTSLIB)

  On Debian/Ubuntu systems you can do this with the command:

  apt-get install git
END


  `which cc` or die <<END;
  'cc' command not in path. Please install it and try again.
  (or to skip Bio::DB::HTS/htslib install re-run with --NO_HTSLIB)

  On Debian/Ubuntu systems you can do this with the command:

  apt-get install build-essential
END

  `which make` or die <<END;
  'make' command not in path. Please install it and try again.
  (or to skip Bio::DB::HTS/htslib install re-run with --NO_HTSLIB)

  On Debian/Ubuntu systems you can do this with the command:

  apt-get install build-essential
END
  
  # List the required libraries with their packages
  my %libs = ( 
    'zlib.h' =>  'zlib1g-dev', 
    'lzma.h' =>  'liblzma-dev', 
    'bzlib.h' => 'libbz2-dev'
  );

  my $msg = '';
  my $this_os =  $^O;

  if ($this_os ne 'darwin' ) {
 
    my $default_msg = qq{%s library header(s) not found in /usr/include. Please install it and try again.
(or to skip Bio::DB::HTS/htslib install re-run with --NO_HTSLIB)

On Debian/Ubuntu systems you can do this with the command:

apt-get install %s};
    my @missing_header = ();
    my @missing_library = (); 
    # Loop over the required libraries
    foreach my $lib (sort(keys(%libs))) {
      unless(-e '/usr/include/'.$lib){
        push(@missing_header, $lib);
	push(@missing_library, $libs{$lib});
      }
    }
    my $header_string = join( ', ', @missing_header);
    my $install_string = join( ' ', @missing_library);
    die(sprintf($default_msg, $header_string, $install_string). "\n\n") if($header_string ne '');
  }
    
  # STEP 1: Create a clean directory for building
  my $htslib_install_dir = $LIB_DIR;
  my $curdir = getcwd;
  chdir $htslib_install_dir;
  my $actualdir = getcwd;

  # STEP 2: Check out HTSLIB / or make this a download?
  print(" - checking out HTSLib\n");
  system "git clone -b $HTSLIB_VERSION https://github.com/samtools/htslib.git";
  -d './htslib' or die "git clone seems to have failed. Could not find $htslib_install_dir/htslib directory";
  chdir './htslib';

  # Step 3: Build libhts.a
  print(" - building HTSLIB in $htslib_install_dir/htslib\n");
  print( "In ".getcwd."\n" );
  # patch makefile
  rename 'Makefile','Makefile.orig' or die "Couldn't rename Makefile to Makefile.orig: $!";
  open my $in, '<','Makefile.orig'     or die "Couldn't open Makefile for reading: $!";
  open my $out,'>','Makefile.new' or die "Couldn't open Makefile.new for writing: $!";

  while (<$in>) {
    chomp;
    if (/^CFLAGS/ && !/-fPIC/) {
      s/#.+//;  # get rid of comments
      $_ .= " -fPIC -Wno-unused -Wno-unused-result";
    }
  }
  continue {
    print $out $_,"\n";
  }

  close $in;
  close $out;
  rename 'Makefile.new','Makefile' or die "Couldn't rename Makefile.new to Makefile: $!";
  system "make";
  -e 'libhts.a' or die "Compile didn't complete. No libhts.a library file found";

  chdir $curdir;
  my $retval = Cwd::realpath("$htslib_install_dir/htslib") ;
}


# INSTALL Bio::DB::HTS
######################
sub install_biodbhts() {

  print "Attempting to install Bio::DB::HTS and htslib.\n\n>>> If this fails, try re-running with --NO_HTSLIB\n\n";

  my $htslib_location = install_htslib();
  rmtree( $DEST_DIR.'/tmp' );

  #Now install Bio::DB::HTS proper
  my $biodbhts_github_url = "https://github.com/Ensembl/Bio-DB-HTS";
  my $biodbhts_zip_github_url = "$biodbhts_github_url/archive/$BIOHTS_VERSION.zip";
  my $biodbhts_zip_download_file = $DEST_DIR.'/tmp/biodbhts.zip';

  mkdir $DEST_DIR unless -d $DEST_DIR;
  mkdir $DEST_DIR.'/tmp';
  download_to_file($biodbhts_zip_github_url, $biodbhts_zip_download_file);
  print " - unpacking $biodbhts_zip_download_file to $DEST_DIR/tmp/\n" unless $QUIET;
  unpack_arch($biodbhts_zip_download_file, "$DEST_DIR/tmp/");

  my $tmp_name = -d "$DEST_DIR/tmp/Bio-HTS-$BIOHTS_VERSION" ? "Bio-HTS-$BIOHTS_VERSION" : "Bio-DB-HTS-$BIOHTS_VERSION";

  print "$DEST_DIR/tmp/$tmp_name - moving files to $BIODBHTS_DIR\n" unless $QUIET;
  rmtree($BIODBHTS_DIR);
  move("$DEST_DIR/tmp/$tmp_name", $BIODBHTS_DIR) or die "ERROR: Could not move directory\n".$!;

  print( " - making Bio::DB:HTS\n" );
  # patch makefile
  chdir $BIODBHTS_DIR;
  system "perl Build.PL --htslib $htslib_location";
  system "./Build";
  chdir ".";

  #move the library
  my $pdir = getcwd;

  #Perl modules to go alongside the API
  dircopy("lib/Bio",$REALPATH_DEST_DIR);

  #The shared object XS library
  if( -e "blib/arch/auto/Bio/DB/HTS/HTS.so" ) {
    copy( "blib/arch/auto/Bio/DB/HTS/Faidx/Faidx.so", "..")
      or die "ERROR: Could not copy shared Faidx.so library:$!\n";
    copy( "blib/arch/auto/Bio/DB/HTS/HTS.so", "..")
      or die "ERROR: Could not copy shared HTS.so library:$!\n";
  }
  elsif( -e "blib/arch/auto/Bio/DB/HTS/HTS.bundle" ) {
    copy( "blib/arch/auto/Bio/DB/HTS/Faidx/Faidx.bundle", "..")
      or die "ERROR: Could not copy shared Faidx.bundle library:$!\n";
    copy( "blib/arch/auto/Bio/DB/HTS/HTS.bundle", "..")
      or die "ERROR: Could not copy shared HTS.bundle library:$!\n";
  }
  else {
    die "ERROR: Shared Bio::DB:HTS library not found\n";
  }

  chdir $pdir;
}

sub dircopy {
  my ($from, $to) = @_;

  opendir FROM, $from;

  foreach my $file(grep {!/^\.\.?/} readdir FROM) {

    # dir?
    if(-d "$from/$file") {
      mkdir("$to/$file") unless -d "$to/$file";
      dircopy("$from/$file", "$to/$file");
    }
    else {
      copy("$from/$file", "$to/$file");
    }
  }

  closedir FROM;
}


# INSTALL BIOPERL
#################
sub bioperl() {

  # now get BioPerl
  print " - fetching BioPerl\n" unless $QUIET;

  my $bioperl_file = (split /\//, $BIOPERL_URL)[-1];

  my $target_file = $DEST_DIR.'/tmp/'.$bioperl_file;

  download_to_file($BIOPERL_URL, $target_file);

  print " - unpacking $target_file\n" unless $QUIET;
  unpack_arch("$DEST_DIR/tmp/$bioperl_file", "$DEST_DIR/tmp/");

  print " - moving files\n" unless $QUIET;

  my $bioperl_dir;

  if($BIOPERL_URL =~ /github/) {
    $bioperl_file =~ s/\.zip//;
    $bioperl_dir = "bioperl-live-".$bioperl_file;
  }
  else {
    $bioperl_file =~ /(bioperl.+?)\.tar\.gz/i;
    $bioperl_dir = $1;
  }

  opendir BIO, "$DEST_DIR/tmp/$bioperl_dir/Bio/";
  move("$DEST_DIR/tmp/$bioperl_dir/Bio/$_", "$DEST_DIR/$_") for readdir BIO;
  closedir BIO;

  rmtree("$DEST_DIR/tmp") or die "ERROR: Failed to remove directory $DEST_DIR/tmp\n";
}


# TEST
######
sub test() {

  print "\nTesting VEP installation\n" unless $QUIET;

  eval q{use Test::Harness; use Test::Exception; };
  if(!$@) {
    opendir TEST, "$dirname\/t";
    my @test_files = map {"$dirname\/t\/".$_} grep {!/^\./ && /\.t$/} readdir TEST;
    closedir TEST;

    # haplo.pl tests require Set::IntervalTree which is not installed here
    eval q{use Set::IntervalTree};
    @test_files = grep {!/Haplo/} @test_files if $@;

    print "Warning: Tests failed, VEP may not run correctly\n" unless runtests(@test_files);
  }
  else {
    my $test_vep = `perl -I $DEST_DIR $dirname/vep --help 2>&1`;

    $test_vep =~ /ENSEMBL VARIANT EFFECT PREDICTOR/ or die "ERROR: Testing VEP script failed with the following error\n$test_vep\n";
  }

  print " - OK!\n" unless $QUIET;
}


# CACHE FILES
#############
sub cache() {

  my $ok;

  if($AUTO) {
    $ok = $AUTO =~ /c/i ? 'y' : 'n';
  }
  else {
    print "\nThe VEP can either connect to remote or local databases, or use local cache files.\nUsing local cache files is the fastest and most efficient way to run the VEP\n" unless $QUIET;
    print "Cache files will be stored in $CACHE_DIR\n" unless $QUIET;

    print "Do you want to install any cache files (y/n)? ";

    $ok = <>;
  }

  if($ok !~ /^y/i) {
    print "Skipping cache installation\n" unless $QUIET;
    return;
  }

  # check cache dir exists
  if(!(-e $CACHE_DIR)) {
    if(!$AUTO) {
      print "Cache directory $CACHE_DIR does not exists - do you want to create it (y/n)? ";

      my $ok = <>;

      if($ok !~ /^y/i) {
        print "Exiting\n";
        exit(0);
      }
    }

    mkdir($CACHE_DIR) or die "ERROR: Could not create directory $CACHE_DIR\n";
  }

  mkdir($CACHE_DIR.'/tmp') unless -e $CACHE_DIR.'/tmp';

  # get list of species
  print " - getting list of available cache files\n" unless $QUIET;

  my $bgzip = `which bgzip`;
  chomp($bgzip);
  $bgzip ||= "$HTSLIB_DIR/bgzip";

  my $tabix = `which tabix`;
  chomp($tabix);
  $tabix ||= "$HTSLIB_DIR/tabix";
  
  my $num = 1;
  my $species_list;
  my $URL_TO_USE = (-e $tabix) ? $CACHE_URL_INDEXED : $CACHE_URL;

  if($URL_TO_USE =~ /^ftp/i) {
    $URL_TO_USE =~ m/(ftp:\/\/)?(.+?)\/(.+)/;
    $ftp = Net::FTP->new($2, Passive => 1) or die "ERROR: Could not connect to FTP host $2\n$@\n";
    $ftp->login($FTP_USER) or die "ERROR: Could not login as $FTP_USER\n$@\n";
    $ftp->binary();

    foreach my $sub(split /\//, $3) {
      $ftp->cwd($sub) or die "ERROR: Could not change directory to $sub\n$@\n";
    }

    push @files, grep {$_ =~ /tar.gz/} $ftp->ls;
  }
  else {
    opendir DIR, $URL_TO_USE;
    @files = grep {$_ =~ /tar.gz/} readdir DIR;
    closedir DIR;
  }

  # if we don't have a species list, we'll have to guess
  if(!scalar(@files)) {
    print "Could not get current species list - using predefined list instead\n";
    print "For more species, see http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#pre\n";

    @files = (
      "bos_taurus_vep_".$DATA_VERSION."_UMD3.1.tar.gz",
      "danio_rerio_vep_".$DATA_VERSION."_GRCz11.tar.gz",
      "homo_sapiens_vep_".$DATA_VERSION."_GRCh37.tar.gz",
      "homo_sapiens_vep_".$DATA_VERSION."_GRCh38.tar.gz",
      "mus_musculus_vep_".$DATA_VERSION."_GRCm38.tar.gz",
      "rattus_norvegicus_vep_".$DATA_VERSION."_Rnor_6.0.tar.gz",
    );
  }
  else {
    # sort
    @files = sort {($a =~ /homo_sapiens/) <=> ($b =~ /homo_sapiens/) || $a cmp $b} @files;
  }

  foreach my $file(@files) {
    $species_list .= $num++." : ".$file."\n";
  }

  if($AUTO) {
    if($SPECIES->[0] eq 'all') {
      @indexes = (1..(scalar @files));
    }

    else {
      foreach my $sp(@$SPECIES) {
        my @matches;

        for my $i(0..$#files) {
          if($sp =~ /refseq|merged/i) {
            push @matches, $i + 1 if $files[$i] =~ /$sp/i;
          }
          else {
            push @matches, $i + 1 if $files[$i] =~ /$sp/i && $files[$i] !~ /refseq|merged/i;
          }
        }

        # grep assembly if supplied
        @matches = grep {$files[$_ - 1] =~ /\_$ASSEMBLY\./} @matches if $ASSEMBLY;

        if(scalar @matches == 1) {
          push @indexes, @matches;
        }
        elsif(scalar @matches > 1) {
          # xenopus_tropicalis_vep_76_JGI_4.2.tar.gz

          my @assemblies = ();
          foreach my $m(@matches) {
            $files[$m-1] =~ m/\_vep\_$DATA_VERSION\_(.+?)\.tar\.gz/;
            push @assemblies, $1 if $1;
          }

          die("ERROR: Multiple assemblies found (".join(", ", @assemblies).") for $sp; select one using --ASSEMBLY [name]\n")
        }
      }
    }

    die("ERROR: No matching species found") unless scalar @indexes;

    # uniquify and sort
    @indexes = sort {$a <=> $b} keys %{{map {$_ => 1} @indexes}};
  }
  else {
    print "The following species/files are available; which do you want (can specify multiple separated by spaces or 0 for all): \n$species_list\n? ";
    @indexes = split /\s+/, <>;

    # user wants all species found
    if(scalar @indexes == 1 && $indexes[0] == 0) {
      @indexes = 1..(scalar @files);
    }
  }

  foreach my $file(@indexes) {
    die("ERROR: File number $file not valid\n") unless defined($file) && $file =~ /^[0-9]+$/ && defined($files[$file - 1]);

    my $file_path = $files[$file - 1];

    my $refseq = 0;
    my ($species, $assembly, $file_name);

    if($file_path =~ /\//) {
      ($species, $file_name) = (split /\//, $file_path);
      $file_name =~ m/^(\w+?)\_vep\_\d+\_(.+?)\.tar\.gz/;
      $assembly = $2;
    }
    else {
      $file_name = $file_path;
      $file_name =~ m/^(\w+?)\_vep\_\d+\_(.+?)\.tar\.gz/;
      $species = $1;
      $assembly = $2;
    }

    push @store_species, $species;

    # check if user already has this species and version
    if(-e "$CACHE_DIR/$species/$DATA_VERSION\_$assembly") {

      my $ok;

      print "\nWARNING: It looks like you already have the cache for $species $assembly (v$DATA_VERSION) installed.\n" unless $QUIET;

      if($AUTO) {
        print "\nDelete the folder $CACHE_DIR/$species/$DATA_VERSION\_$assembly and re-run INSTALL.pl if you want to re-install\n";
      }
      else {
        print "If you continue the existing cache will be overwritten.\nAre you sure you want to continue (y/n)? ";

        $ok = <>;
      }

      if($ok !~ /^y/i) {
        print " - skipping $species\n" unless $QUIET;
        next;
      }

      rmtree("$CACHE_DIR/$species/$DATA_VERSION\_$assembly") or die "ERROR: Could not delete directory $CACHE_DIR/$species/$DATA_VERSION\_$assembly\n";
    }

    if($species =~ /refseq/i) {
      print "NB: Remember to use --refseq when running the VEP with this cache!\n" unless $QUIET;
    }
    if($species =~ /merged/i) {
      print "NB: Remember to use --merged when running the VEP with this cache!\n" unless $QUIET;
    }

    my $target_file = "$CACHE_DIR/tmp/$file_name";
    if($URL_TO_USE =~ /^ftp/) {
      print " - downloading $URL_TO_USE/$file_path\n" unless $QUIET;
      if(!$TEST) {
        $ftp->get($file_name, $target_file) or download_to_file("$URL_TO_USE/$file_path", $target_file);

        my $checksums = "CHECKSUMS";
        my $checksums_target_file = "$CACHE_DIR/tmp/$checksums";
        $ftp->get($checksums, $checksums_target_file) or download_to_file("$URL_TO_USE/$checksums", $checksums_target_file);
        if (-e $checksums_target_file) {
          my $sum_download = `sum $target_file`;
          $sum_download =~ m/([0-9]+)(\s+)([0-9]+)/;
          my $checksum_download = $1;
          $checksum_download =~ s/^0*//;
          my $sum_ftp = `grep $file_name $checksums_target_file`;
          $sum_ftp =~ s/^0*//;
          if ($sum_download && $sum_ftp) {
            die("ERROR: checksum for $target_file doesn't match checksum in CHECKSUMS file on FTP site\n") if ($sum_ftp !~ m/^$checksum_download\s+/);
          }
        }
      }
    }
    else {
      print " - copying $URL_TO_USE/$file_path\n" unless $QUIET;
      copy("$URL_TO_USE/$file_path", $target_file) unless $TEST;
    }

    print " - unpacking $file_name\n" unless $QUIET;


    unpack_arch($target_file, $CACHE_DIR.'/tmp/') unless $TEST;

    # does species dir exist?
    if(!-e "$CACHE_DIR/$species" && !$TEST) {
      mkdir("$CACHE_DIR/$species") or die "ERROR: Could not create directory $CACHE_DIR/$species\n";
    }

    # move files
    unless($TEST) {
      opendir CACHEDIR, "$CACHE_DIR/tmp/$species/";
      move("$CACHE_DIR/tmp/$species/$_", "$CACHE_DIR/$species/$_") for readdir CACHEDIR;
      closedir CACHEDIR;
    }
    
    if(((-e $bgzip && -e $tabix) || $CONVERT) && !$TEST) {
      unless($QUIET) {
        print " - converting cache, this may take some time but will allow VEP to look up variants and frequency data much faster\n";
        print " - use CTRL-C to cancel if you do not wish to convert this cache now (you may run convert_cache.pl later)\n";
      }
      system("perl $dirname/convert_cache.pl --dir $CACHE_DIR --species $species --version $DATA_VERSION\_$assembly --bgzip $bgzip --tabix $tabix") == 0 or print STDERR "WARNING: Failed to run convert script\n";
    }
  }
}


# FASTA FILES
#############
sub fasta() {

  ### SPECIAL CASE GRCh37
  if((grep {$files[$_ - 1] =~ /GRCh37/} @indexes) || (defined($ASSEMBLY) && $ASSEMBLY eq 'GRCh37')) {

    # can't install other species at same time as the FASTA URL has to be changed
    if(grep {$files[$_ - 1] !~ /GRCh37/} @indexes) {
      die("ERROR: For technical reasons this installer is unable to install GRCh37 caches alongside others; please install them separately\n");
    }

    # change URL to point to last e! version that had GRCh37 downloads
    elsif($FASTA_URL =~ /ftp/) {
      print "\nWARNING: Changing FTP URL for GRCh37\n";
      $FASTA_URL =~ s/$DATA_VERSION/75/;
    }
  }

  my $ok;

  if($AUTO) {
    $ok = $AUTO =~ /f/i ? 'y' : 'n';
  }
  else {
    print "\nThe VEP can use FASTA files to retrieve sequence data for HGVS notations and reference sequence checks.\n" unless $QUIET;
    print "FASTA files will be stored in $CACHE_DIR\n" unless $QUIET;
    print "Do you want to install any FASTA files (y/n)? ";

    $ok = <>;
  }

  if($ok !~ /^y/i) {
    print "Skipping FASTA installation - Exiting\n";
    return;
  }

  my @dirs = ();

  if($FASTA_URL =~ /^ftp/i) {
    $FASTA_URL =~ m/(ftp:\/\/)?(.+?)\/(.+)/;
    $ftp = Net::FTP->new($2, Passive => 1) or die "ERROR: Could not connect to FTP host $2\n$@\n";
    $ftp->login($FTP_USER) or die "ERROR: Could not login as $FTP_USER\n$@\n";
    $ftp->binary();

    foreach my $sub(split /\//, $3) {
      $ftp->cwd($sub) or die "ERROR: Could not change directory to $sub\n$@\n";
    }

    push @dirs, sort $ftp->ls;
  }
  else {
    opendir DIR, $FASTA_URL;
    @dirs = grep {-d $FASTA_URL.'/'.$_ && $_ !~ /^\./} readdir DIR;
    closedir DIR;
  }

  my $species_list = '';
  my $num = 1;
  foreach my $dir(@dirs) {
    $species_list .= $num++." : ".$dir."\n";
  }

  my @species;
  if($AUTO) {
    if($SPECIES->[0] eq 'all') {
      @species = scalar @store_species ? @store_species : @dirs;
    }
    else {
      @species = scalar @store_species ? @store_species : @$SPECIES;
    }
  }
  else {
    print "FASTA files for the following species are available; which do you want (can specify multiple separated by spaces, \"0\" to install for species specified for cache download): \n$species_list\n? ";

    my $input = <>;
    my @nums = split /\s+/, $input;

    @species = @store_species if grep {$_ eq '0'} @nums;
    push @species, $dirs[$_ - 1] for grep {$_ > 0} @nums;
  }

  foreach my $species(@species) {

    # remove refseq name
    my $orig_species = $species;
    $species =~ s/_refseq//;
    $species =~ s/_merged//;

    my @files;
    my $dna_path;

    if($ftp) {
      $ftp->cwd($species) or die "ERROR: Could not change directory to $species\n$@\n";
      $ftp->cwd('dna_index') or $ftp->cwd('dna') or die "ERROR: Could not change directory to dna\n$@\n";
      @files = $ftp->ls;
      $dna_path = $ftp->pwd =~ /dna_index/ ? 'dna_index' : 'dna';
    }
    else {
      if(!opendir DIR, "$FASTA_URL/$species/dna") {
        warn "WARNING: Could not read from directory $FASTA_URL/$species/dna\n$@\n";
        next;
      }
      @files = grep {$_ !~ /^\./} readdir DIR;
      closedir DIR;
    }

    # remove repeat/soft-masked files
    @files = grep {!/_(s|r)m\./} @files;

    my ($file) = grep {/primary_assembly.fa.gz$/} @files;
    ($file) = grep {/toplevel.fa.gz$/} @files if !defined($file);

    unless(defined($file)) {
      warn "WARNING: No download found for $species\n";
      next;
    }

    # work out assembly version from file name
    my $uc_species = ucfirst($species);
    $file =~ m/^$uc_species\.(.+?)(\.)?(\d+)?\.dna/;
    my $assembly = $1;

    # second number could be an Ensembl release number (pre-76) or part of the assembly name
    if(defined($3)) {
      if(!grep {$3 == $_} (69..75)) {
        $assembly .= $2.$3;
      }
    }

    die("ERROR: Unable to parse assembly name from $file\n") unless $assembly;

    my $ex = "$CACHE_DIR/$orig_species/$DATA_VERSION\_$assembly/$file";
    my $ex_unpacked = $ex;
    $ex_unpacked =~ s/\.gz$//;
    if(-e $ex || -e $ex_unpacked) {
      print "Looks like you already have the FASTA file for $orig_species, skipping\n" unless $QUIET;

      if($ftp) {
        $ftp->cwd('../');
        $ftp->cwd('../');
      }
      next;
    }

    # create path
    mkdir($CACHE_DIR) unless -d $CACHE_DIR || $TEST;
    mkdir("$CACHE_DIR/$orig_species") unless -d "$CACHE_DIR/$orig_species" || $TEST;
    mkdir("$CACHE_DIR/$orig_species/$DATA_VERSION\_$assembly") unless -d "$CACHE_DIR/$orig_species/$DATA_VERSION\_$assembly" || $TEST;

    if($ftp) {
      print " - downloading $file\n" unless $QUIET;
      if(!$TEST) {
        $ftp->get($file, $ex) or download_to_file("$FASTA_URL/$species/$dna_path/$file", $ex);
      }
    }
    else {
      print " - copying $file\n" unless $QUIET;
      copy("$FASTA_URL/$species/dna/$file", $ex) unless $TEST;
    }

    my $bgzip = `which bgzip`;
    chomp($bgzip);
    $bgzip ||= "$HTSLIB_DIR/bgzip";

    eval q{ use Bio::DB::HTS::Faidx; };
    my $can_use_faidx = $@ ? 0 : 1;

    if($can_use_faidx) {

      if($dna_path !~ /dna_index/ && -e $bgzip && $CAN_USE_GZIP) {
        print " - converting sequence data to bgzip format, this may take some time...\n" unless $QUIET;
        my $curdir = getcwd;
        my $bgzip_convert = "gzip -dc $ex | $bgzip -c > $ex\.bgz; mv $ex\.bgz $ex";
        my $bgzip_result = `$bgzip_convert` unless $TEST;

        if( $? != 0 ) {
          die "FASTA gzip to bgzip conversion failed: $bgzip_result\n" unless $TEST;
        }
        else {
          print " - conversion successful\n";
        }
      }

      my $got_indexes = 0;
      foreach my $index_file(grep {/$file\..+/} @files) {
        if(!$TEST) {
          $index_file =~ /$file(\..+)/;
          print " - downloading $index_file\n" unless $QUIET;
          $ftp->get($index_file, $ex.$1) or download_to_file("$FASTA_URL/$species/$dna_path/$index_file", $ex.$1);
          $got_indexes++;
        }
      }

      unless($got_indexes == 2) {
        print " - indexing FASTA file\n" unless $QUIET;
        Bio::DB::HTS::Faidx->new($ex) unless $TEST;
        print " - indexing OK\n" unless $QUIET;
      }

      print "\nThe FASTA file should be automatically detected by the VEP when using --cache or --offline.\nIf it is not, use \"--fasta $ex\"\n\n" unless $QUIET;
    }

    elsif($NO_HTSLIB) {
      print " - extracting data\n" unless $QUIET;
      unpack_arch($ex, "$CACHE_DIR/$orig_species/$DATA_VERSION\_$assembly/") unless $TEST;

      print " - attempting to index\n" unless $QUIET;
      eval q{
        use Bio::DB::Fasta;
      };
      if($@) {
        print "Indexing failed - VEP will attempt to index the file the first time you use it\n" unless $QUIET;
      }
      else {
        print " - indexing FASTA file\n" unless $QUIET;
        Bio::DB::Fasta->new($ex_unpacked) unless $TEST;
        print " - indexing OK\n" unless $QUIET;
      }

      print "\nThe FASTA file should be automatically detected by the VEP when using --cache or --offline.\nIf it is not, use \"--fasta $ex_unpacked\"\n\n" unless $QUIET;
    }

    if($ftp) {
      $ftp->cwd('../');
      $ftp->cwd('../');
    }
  }
}

# PLUGINS
#########

sub plugins() {
  my $ok;

  if($AUTO) {
    $ok = $AUTO =~ /p/i ? 'y' : 'n';
  }
  else {
    print "\nThe VEP can use plugins to add functionality and data.\n" unless $QUIET;
    print "Plugins will be installed in $PLUGINS_DIR\n" unless $QUIET;

    print "Do you want to install any plugins (y/n)? ";

    $ok = <>;
  }

  if($ok !~ /^y/i) {
    print "Skipping plugin installation\n" unless $QUIET;
    return;
  }

  # check plugin installation dir exists
  if(!(-e $PLUGINS_DIR)) {
    if(!$AUTO) {
      print "Plugins directory $PLUGINS_DIR does not exists - do you want to create it (y/n)? ";

      my $ok = <>;

      if($ok !~ /^y/i) {
        print "Exiting\n";
        exit(0);
      }
    }

    mkpath($PLUGINS_DIR) or die "ERROR: Could not create directory $PLUGINS_DIR\n";
  }
  mkpath($PLUGINS_DIR) unless -e $PLUGINS_DIR;

  my $plugin_url_root = $PLUGIN_URL.'/release/'.$API_VERSION;

  # download and eval plugin config file
  my $plugin_config_file = $PLUGINS_DIR.'/plugin_config.txt';
  download_to_file($plugin_url_root.'/plugin_config.txt', $plugin_config_file);

  die("ERROR: Could not access plugin config file $plugin_config_file\n") unless($plugin_config_file && -e $plugin_config_file);

  open IN, $plugin_config_file;
  my @content = <IN>;
  close IN;

  my $VEP_PLUGIN_CONFIG = eval join('', @content);
  die("ERROR: Could not eval VEP plugin config file: $@\n") if $@;
  my @plugins = @{$VEP_PLUGIN_CONFIG->{plugins}};

  # get sections
  my @sections = grep {defined($_)} map {defined($_->{section}) ? $_->{section} : undef} @plugins;

  # unique sort in same order
  my ($prev, @new);
  for(@sections) {
    if(!defined($prev) || $prev ne $_) {
      push @new, $_;
    }
    $prev = $_;
  }
  @sections = @new;
  push @sections, '';

  # generate list to present to user
  my (%by_number, %by_key);
  my $i = 1;
  my $plugin_list = '';

  # get length of longest label/key
  my $length = length((sort {length($a->{key}) <=> length($b->{key})} @plugins)[-1]->{key});

  foreach my $section(@sections) {
    my $section_name = $section || (scalar @sections > 1 ? 'Other plugins' : 'Plugins');

    my @section_plugins;

    # check that plugins have plugin_url defined
    # otherwise we can't download it
    if($section eq '') {
      @section_plugins = grep {$_->{plugin_url} && !defined($_->{section})} @plugins;
    }
    else {
      @section_plugins = grep {$_->{plugin_url} && defined($_->{section}) && $_->{section} eq $section} @plugins;
    }

    next unless scalar @section_plugins;
    $plugin_list .= "# $section_name\n";

    $_->{plugin_number} = $i++ for @section_plugins;
    $by_number{$_->{plugin_number}} = $_ for @section_plugins;
    $by_key{lc($_->{key})} = $_ for @section_plugins;

    $plugin_list .= sprintf(
      "%4i: %*s - %s\n",
      $_->{plugin_number},
      0 - $length,
      $_->{key},
      $_->{helptip} || ''
    ) for @section_plugins;
  }

  # now establish which we are installing
  my (@indexes, @selected_plugins);

  # either from user input
  if(!$AUTO) {
    print "\nThe following plugins are available; which do you want (can specify multiple separated by spaces or 0 for all): \n$plugin_list\n? ";
    @indexes = split /\s+/, <>;

    # user wants all species found
    if(scalar @indexes == 1 && $indexes[0] == 0) {
      @indexes = 1..(scalar keys %by_number);
    }

    @selected_plugins = map {$by_number{$_}} grep {$by_number{$_}} @indexes;
  }

  # or from list passed on command line
  else {
    if(lc($PLUGINS->[0]) eq 'all' || $PLUGINS->[0] eq '0') {
      @selected_plugins = sort {$a->{key} cmp $b->{key}} values %by_key;
    }
    else {
      @selected_plugins = map {$by_key{lc($_)}} grep {$by_key{lc($_)}} @$PLUGINS;
    }

    my @not_found = grep {!$by_key{lc($_)}} @$PLUGINS;
    if(@not_found) {
      printf(
        "\nWARNING: The following plugins have not been found: %s\nAvailable plugins: %s\n",
        join(",", @not_found),
        join(",", sort map {$_->{key}} values %by_key)
      );
    }

    if(!@selected_plugins) {
      printf("\nERROR: No valid plugins given\n");
      return;
    }
  }

  # store a flag to warn user at end if any plugins require additional setup
  my $requires_install_or_data = 0;

  foreach my $pl(@selected_plugins) {
    printf("\n - installing \"%s\"\n", $pl->{key});

    my $local_file = $PLUGINS_DIR.'/'.$pl->{key}.'.pm';

    # overwrite?
    if(-e $local_file) {
      printf(
        "%s already installed; %s",
        $pl->{key},
        $AUTO ? "overwriting\n" : "do you want to overwrite (probably OK if updating) (y/n)? "
      );

      my $ok = $AUTO ? 'y' : <>;

      if($ok !~ /^y/i) {
        print " - Skipping\n";
        next;
      }
    }

    # download
    download_to_file($pl->{plugin_url}, $local_file);

    # warn if failed
    unless(-e $local_file) {
      print " - WARNING: Failed to download/install ".$pl->{key}."\n";
      next;
    }

    # additional setup required?
    if($pl->{requires_install} || $pl->{requires_data}) {
      print " - This plugin requires installation\n" if $pl->{requires_install};
      print " - This plugin requires data\n" if $pl->{requires_data};
      print " - See $local_file for details\n";

      $requires_install_or_data = 1;
    }

    else {
      printf(
        " - add \"--plugin %s%s\" to your VEP command to use this plugin\n",
        $pl->{key},
        $pl->{params} ? ',[options]' : ''
      );
    }

    print " - OK\n";
  }

  print "\nNB: One or more plugins that you have installed will not work without installation or downloading data; see logs above\n" if $requires_install_or_data;
}

# OTHER SUBS
############

sub download_to_file {
  my ($url, $file) = @_;

  $url =~ s/([a-z])\//$1\:21\// if $url =~ /ftp/ && $url !~ /\:21/;

  if($CAN_USE_CURL) {
    my $response = `curl -s -o $file -w '%{http_code}' --location "$url" `;
    if ( $response != 200 && $response != 226) {
      print "curl failed ($response), trying to fetch using LWP::Simple\n" unless $QUIET;
      $CAN_USE_CURL = 0;
      download_to_file($url, $file);
    }
  }

  elsif($CAN_USE_LWP) {
    my $response = getstore($url, $file);

    unless($response == 200) {

      # try no proxy
      $ua->no_proxy('github.com');

      $response = getstore($url, $file);

      unless($response == 200) {
        print "LWP::Simple failed ($response), trying to fetch using HTTP::Tiny\n" unless $QUIET;
        $CAN_USE_LWP = 0;
        download_to_file($url, $file);
      }
    }
  }
  elsif($CAN_USE_HTTP_TINY) {
    my $response = HTTP::Tiny->new(no_proxy => 'github.com')->get($url);

    if($response->{success}) {
      open OUT, ">$file" or die "Could not write to file $file\n";
      binmode OUT;
      print OUT $response->{content};
      close OUT;
    }
    else {
      #warn "WARNING: Failed to fetch from $url\nError code: $response->{reason}\nError content:\n$response->{content}\nTrying without no_proxy\n" unless $QUIET;
      $response = HTTP::Tiny->new->get($url);

      if($response->{success}) {
        open OUT, ">$file" or die "Could not write to file $file\n";
        binmode OUT;
        print OUT $response->{content};
        close OUT;
      }
      else {
        die("ERROR: Failed last resort of using HTTP::Tiny to download $url\n");
      }
    }
  }
  else {
    die("ERROR: Unable to download files without curl, LWP or HTTP::Tiny installed\n");
  }
}

# unpack a tarball
sub unpack_arch {
  my ($arch_file, $dir) = @_;

  if($CAN_USE_ARCHIVE) {
    my $ar = Archive::Extract->new(archive => $arch_file);
    my $ok = $ar->extract(to => $dir) or die $ar->error;
  }
  else {
    if($arch_file =~ /.zip$/ && $CAN_USE_UNZIP) {
      `unzip -qq $arch_file -d $dir` and die("ERROR: Failed to unpack file $arch_file\n");
    }
    elsif($arch_file =~ /(\.tar\.|\.t)gz$/ && $CAN_USE_TAR) {
      `tar -C $dir -zxf $arch_file` and die("ERROR: Failed to unpack file $arch_file\n");
    }
    elsif($arch_file =~ /\.gz$/ && $arch_file !~ /(\.tar\.|\.t)gz$/ && $CAN_USE_GZIP) {
      my $unpacked = $arch_file;
      $unpacked =~ s/.*\///g;
      $unpacked =~ s/\.gz$//;
      `gzip -dc $arch_file > $dir/$unpacked` and die("ERROR: Failed to unpack file $arch_file\n");
    }
    else {
      die("ERROR: Unable to unpack file $arch_file without Archive::Extract or tar/unzip/gzip\n");
    }
  }

  unlink($arch_file);
}

sub usage {
  my $versions = get_version_string();

  my $usage =<<END;
#---------------#
# VEP INSTALLER #
#---------------#

versions
  $versions

http://www.ensembl.org/info/docs/variation/vep/vep_script.html#installer

Usage:
perl INSTALL.pl [arguments]

Options
=======

-h | --help        Display this message and quit

-d | --DESTDIR     Set destination directory for API install (default = './')
--CACHE_VERSION    Set data (cache, FASTA) version to install if different from --VERSION (default = $VERSION)
-c | --CACHEDIR    Set destination directory for cache files (default = '$ENV{HOME}/.vep/')

-a | --AUTO        Run installer without user prompts. Use "a" (API + Faidx/htslib),
                   "l" (Faidx/htslib only), "c" (cache), "f" (FASTA), "p" (plugins) to specify
                   parts to install e.g. -a ac for API and cache
-n | --NO_UPDATE   Do not check for updates to ensembl-vep
-s | --SPECIES     Comma-separated list of species to install when using --AUTO
-y | --ASSEMBLY    Assembly name to use if more than one during --AUTO
-g | --PLUGINS     Comma-separated list of plugins to install when using --AUTO
-r | --PLUGINSDIR  Set destination directory for VEP plugins files (default = '$ENV{HOME}/.vep/Plugins/')
-q | --QUIET       Don't write any status output when using --AUTO
-p | --PREFER_BIN  Use this if the installer fails with out of memory errors
-l | --NO_HTSLIB   Don't attempt to install Faidx/htslib
--NO_BIOPERL       Don't install BioPerl

-t | --CONVERT     Convert downloaded caches to use tabix for retrieving
                   co-located variants (requires tabix)


-u | --CACHEURL    Override default cache URL; this may be a local directory or
                   a remote (e.g. FTP) address.
-f | --FASTAURL    Override default FASTA URL; this may be a local directory or
                   a remote (e.g. FTP) address. The FASTA URL/directory must have
                   gzipped FASTA files under the following structure:
                   [species]/[dna]/
END

  print $usage;
}
