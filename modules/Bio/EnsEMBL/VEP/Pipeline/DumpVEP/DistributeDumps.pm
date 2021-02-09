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
package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::DistributeDumps;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::VEP::Pipeline::DumpVEP::BaseVEP);

use FileHandle;
use File::Path qw(make_path);
use File::Spec;

sub run {
  my $self = shift;
  my $version = $self->param('eg_version') || $self->required_param('ensembl_release');
  my $dir;
  my $division = $self->param('division') || [];
  my @division = ( ref($division) eq 'ARRAY' ) ? @$division : ($division);
  # If division param is defined.
  if ( scalar(@division) ) {
    foreach my $division (@division) {
      $dir=$self->required_param('pipeline_dir')."/".$division.'/dumps';
      $self->DistributeProduction($dir);
      $self->DistributeWeb($dir);
    }
  }
  else {
    die "No division found, cannot distribute dumps";
  }
}



#Distribute the Production dumps for the FTP site
sub DistributeProduction {
  my ($self,$dir) = @_;
  make_path("$dir/production") if (!-d "$dir/production");
  opendir (my $dh, $dir) or die $!;
  while (my $content = readdir($dh)) {
    if ($content =~ /collection$/)
    {
      make_path("$dir/production/$content") if (!-d "$dir/production/$content");
      opendir (my $dh_collection, "$dir/$content") or die $!;
      while (my $file_collection = readdir($dh_collection)) {
        if ($file_collection =~ /gz$/ && $file_collection !~ /tabix/) {
          $self->link_file("$dir/$content/$file_collection", "$dir/production/$content/$file_collection");
        }
      }
      $dh_collection->close();
      $self->compute_checksums("$dir/production/$content/");
    }
    elsif ($content =~ /gz$/ && $content !~ /tabix/) {
      $self->link_file("$dir/$content", "$dir/production/$content");
    }
  }
  $dh->close();
  $self->compute_checksums("$dir/production/");
}


# Distribute the Web dumps for the VEP browser tool.
# These dumps are tabix-converted (they contain files that require tabix to be installed to access)
# only species with variant data are tabix-converted
# you get a big speed increase when looking up co-located variants on the web interface using a tabix-based cache vs not
# VEP installer does include an attempt to install Bio::DB::HTS which allows use of the converted dumps but it doesnt have a 100% success rate. We are experiencing issues with OSX ATM, for example
sub DistributeWeb{
  my ($self,$dir) = @_;

  my %copied; ## save list of tabix'ed tar balls copied

  make_path("$dir/web") if (!-d "$dir/web");
  opendir (my $dh, $dir) or die $!;
  while (my $content = readdir($dh)) {
    if ($content =~ /collection$/)
    {
      make_path("$dir/web/$content") if (!-d "$dir/web/$content");
      opendir (my $dh_collection, "$dir/$content") or die $!;
      while (my $file_collection = readdir($dh_collection)) {
        ## only copy non- tabixed set if tabixed set not already copied
        if ($file_collection =~ /gz$/ && $file_collection !~ /tabix/ && ! $copied{$file_collection}) {
          $self->link_file("$dir/$content/$file_collection", "$dir/web/$content/$file_collection");
        }
        elsif ($file_collection =~ /tabix/) {
          my $file_collection_no_tabix = $file_collection;
          $file_collection_no_tabix =~ s/\_tabixconverted//;
          $copied{$file_collection_no_tabix} = 1;
          $self->link_file("$dir/$content/$file_collection", "$dir/web/$content/$file_collection_no_tabix");
        }
      }
      $dh_collection->close();
    }
    elsif ($content =~ /gz$/ && $content !~ /tabix/ && ! $copied{$content}) {
      $self->link_file("$dir/$content", "$dir/web/$content");
    }
    elsif ($content =~ /tabix/) {
      my $file_no_tabix = $content;
      $file_no_tabix =~ s/\_tabixconverted//;
      $copied{$file_no_tabix} = 1;
      $self->link_file("$dir/$content", "$dir/web/$file_no_tabix");
    }
  }
  $dh->close();
}

sub compute_checksums {
  my ($self,$dir) = @_;
  opendir(my $dh, $dir) or die $!;
  my @files = sort {$a cmp $b} readdir($dh);
  closedir($dh) or die $!;
  my @checksums = ();
  foreach my $file (@files) {
    next if $file =~ /^\./;
    next if $file =~ /^CHECKSUM/;
    my $path = File::Spec->catfile($dir, $file);
    next if (-d $path);
    my $checksum = checksum($path);
    push(@checksums, [$checksum, $file]);
  }
  my $fh = FileHandle->new("$dir/CHECKSUMS", 'w');
  foreach my $entry (@checksums) {
    my $line = join("\t", @{$entry});
    print $fh $line, "\n";
  }
  $fh->close();
}

sub checksum {
  my $path = shift;
  my $checksum = `sum $path`;
  die('sum ' . $path . ' failed: ' . $?) if $?;
  $checksum =~ s/\s* $path//xms;
  chomp($checksum);
  return $checksum;
}

1;
