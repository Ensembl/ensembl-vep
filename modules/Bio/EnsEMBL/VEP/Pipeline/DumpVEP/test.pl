use File::Spec;
use FileHandle;
compute_checksums(1,'/hps/nobackup2/production/ensembl/ensprod/release_dumps/release-45/dump_vep_45/protists/dumps/production/protists_alveolata1_collection');



sub compute_checksums {
  my ($self,$dir) = @_;
  $DB::single = 1;
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
