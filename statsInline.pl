#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage = <<END
Usage: statsInline.pl
       statsInline.pl -dir .

After running cleanInline.pl, reports how many reads passed each phase
  of the pipeline, and how many of the parsed ASVs are singletons
  within their samples.
END
  ;

my $dir = ".";
die $usage
  unless GetOptions('dir=s' => \$dir)
  && @ARGV == 0;
die "Not a directory: $dir\n" unless -d $dir;

my @pearLogs = glob("$dir/*.pear.log");
my @filtered = glob("$dir/*.filtered");
my @parsed = glob("$dir/*.parse.tab");

die "No *.pear.log files in $dir\n" if @pearLogs == 0;
die "No *.filtered files in $dir\n" if @filtered == 0;
die "No *.parse.tab files in $dir\n" if @parsed == 0;

my $nReads = 0;
my $nAssembled = 0;
foreach my $file (@pearLogs) {
  open(my $fh, "<", $file) || die "Cannot read $file\n";
  while(my $line = <$fh>) {
    if ($line =~ m!^Assembled reads .*: *([0-9,]+) */ *([0-9,]+) *!) {
      my ($nAssembledThis, $nReadsThis) = ($1,$2);
      $nAssembledThis =~ s/,//g;
      $nReadsThis =~ s/,//g;
      $nAssembled += $nAssembledThis;
      $nReads += $nReadsThis;
    }
  }
  close($fh) || die "Error reading $file\n";
}
print STDERR "Scanned " . scalar(@pearLogs) . " pear logs\n";
die "No reads\n" if $nReads < 1;

my $nQuality = 0;
foreach my $file (@filtered) {
  my $nQualityThis = `grep -c ">" $file`;
  $nQuality += $nQualityThis;
}
print STDERR "Scanned " . scalar(@filtered) . " filtered files\n";

my $nUsed = 0;
my $nSingleton = 0;
foreach my $file (@parsed) {
  open(my $fh, "<", $file) || die "Cannot read $file\n";
  while (my $line = <$fh>) {
    chomp $line;
    my (undef, $n) = split /\t/, $line;
    next if $n eq "count";
    die "Invalid line\n$line\nin $file" unless $n =~ m/^\d+$/;
    $nUsed += $n;
    $nSingleton++ if $n == 1;
  }
  close($fh) || die "Error reading $file\n";
}

print STDERR sprintf(join("\n",
                          "Total reads:  %6.1f M %5.1f%%",
                          "Assembled:    %6.1f M %5.1f%%",
                          "High-quality: %6.1f M %5.1f%%",
                          "Demultiplexed:%6.1f M %5.1f%%",
                          "Singletons:   %6.1f M %5.1f%%",
                          "Non-unique:   %6.1f M %5.1f%%",
                          ""),
                     $nReads / 1e6, 100,
                     $nAssembled / 1e6, 100 * $nAssembled / $nReads,
                     $nQuality / 1e6, 100 * $nQuality / $nReads,
                     $nUsed / 1e6, 100 * $nUsed / $nReads,
                     $nSingleton / 1e6, 100 * $nSingleton / $nReads,
                     ($nUsed - $nSingleton) / 1e6, 100 * ($nUsed - $nSingleton) / $nReads);
