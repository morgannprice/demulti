#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw{sum};
use Digest::MD5 qw{md5_hex};
use FindBin qw{$RealBin};
use lib $RealBin;
use dmUtils qw{ReadTable ReadFastaEntry};
sub compareEsvNames($$);

my $usearch = "usearch";
my $minN = 4;
my $namePre = "";
my $debug;
my $minLen = 200;

my $usage = <<END
Usage: cleanInline.pl -in *.V4.tab -out outprefix
Input files should be from parseInline.pl

For each sample, filters out rare sequences, uses unoise3 from usearch
to remove uncommon noisy sequences and chimeras.  Writes out.tsv, with
a unified table of counts across all samples (with arbitrary ESV names
based on MD5 hashes), and out.fna, with the actual sequence of each
ESV.

Optional arguments:
-zotu fastafile -- read in existing names of zotus from the file
-minN $minN -- ignore sequences with fewer than this reads in a sample
-usearch $usearch -- the usearch executable
-name $namePre -- the prefix for each ESV name
   By default, uses md5 hashes, not prefix number
-primers 4,6,7 -- ignore primers whose names do not end with these
   numbers
-debug -- voluminous output from usearch
-minLen $minLen -- ignore sequences of below this length
END
;

my @infiles;
my $outPre;
my $oldFasta;
my $primerSpec;
die $usage
  unless GetOptions('usearch=s' => \$usearch,
                    'minN=i' => \$minN, 
                    'in=s{,}' => \@infiles,
                    'primers=s' => \$primerSpec,
                    'name=s' => \$namePre,
                    'out=s' => \$outPre,
                    'minLen=i' => \$minLen,
                    'zotu=s' => \$oldFasta,
                    'debug' => \$debug)
  && defined $outPre
  && @ARGV==0;
die "No input files\n$usage" unless @infiles > 0;
die "minN must be at least 1\n" unless $minN >= 1;
if (-x "$RealBin/$usearch") {
  $usearch = "$RealBin/$usearch";
} elsif (-x $usearch) {
  ;
} else {
  die "Cannot find usearch executable $usearch\n";
}

my %primers = (); # primer => 1 if the end of a known primer name
if (defined $primerSpec) {
  die "Invalid primer specifier $primerSpec\n"
    unless $primerSpec =~ m/^[0-9,]+$/;
  my @primers = split /,/, $primerSpec;
  foreach my $primer (@primers) {
    die "Invalid primer specifier $primerSpec\n"
      unless $primer =~ m/^\d+$/ && $primer > 0;
  }
  %primers = map { $_ => 1 } @primers;
}

foreach my $infile (@infiles) {
  die "No such input file: $infile\n" unless -e $infile;
}

my $tmpPre = "/tmp/cleanInline.$$";
my $tmpFna = "$tmpPre.fna";
my $tmpU = "$tmpPre.u";


my %seqToEsv = (); # sequence to ESV name
my %esvToSeq = (); # ESV name to sequence
my %counts = (); # primer_name => index => ESV name => count
my %indexSeen = ();

if (defined $oldFasta) {
  open(my $fh, "<", $oldFasta) || die "Cannot read $oldFasta\n";
  my $state = {};
  while (my ($esv, $seq) = ReadFastaEntry($fh, $state)) {
    if (exists $seqToEsv{$seq}) {
      print STDERR "Warning, sequence for $esv is a duplicate, ignored\n";
    } else {
      $esv =~ s/[ \t].*//;
      die "Duplicate ESV name $esv in $oldFasta"
        if exists $esvToSeq{$esv};
      $seqToEsv{$seq} = $esv;
      $esvToSeq{$esv} = $seq;
    }
  }
  close($fh) || die "Error reading $oldFasta";
  print STDERR "Read " . scalar(keys %esvToSeq) . " ESV names from $oldFasta\n";
}


my $totReadsIn = 0;
my $totReadsOut = 0;

foreach my $infile (@infiles) {
  my $name = $infile; $name = $1 if $infile =~ m!/([^/]+)$!;
  die "Cannot parse index from $infile"
    unless $name =~ m/^([a-zA-Z]+\d+)[._]/
      || $name =~ m!^/([a-zA-Z]+\d+)[.][^/]+$!;
  my $index = $1;
  die "Duplicate index $index from file $infile\n" if exists $indexSeen{$index};
  $indexSeen{$index} = 1;
  my $base = $infile;
  $base =~ s!^.*/!!;
  my @rows = ReadTable($infile, ["primer_name","count","sequence"]);
  @rows = grep { $_->{count} >= $minN && length($_->{sequence}) >= $minLen } @rows;
  if (keys %primers > 0) {
    @rows = grep { my $name = $_->{primer_name};
                   $name =~ m/[a-zA-Z](\d+)$/ || die "Invalid primer_name $name";
                   my $i = $1;
                   exists $primers{$i} } @rows;
  }
  if (@rows == 0) {
    print STDERR "Warning: no rows in $infile with count >= $minN and length >= $minLen and expected primers\n";
    next;
  }
  my $nReadsThis = sum(map $_->{count}, @rows);
  $totReadsIn += $nReadsThis;
  print STDERR "Considering " . scalar(@rows) . " sequences ($nReadsThis reads) in $infile\n";

  # For each primer, handle each sample separately
  my %byprimer = ();
  foreach my $row (@rows) {
    push @{ $byprimer{$row->{primer_name}} }, $row;
  }

  # Within each primer, they should already be sorted by n
  while (my ($primer, $prows) = each %byprimer) {
    open(my $fhFna, ">", $tmpFna) || die "Cannot write to $tmpFna\n";
    my %numToRow = (); # numbering of the sequences given to usearch
    my $nSeq = 0;
    foreach my $row (@$prows) {
      $numToRow{$nSeq} = $row;
      print $fhFna ">SEQ" . $nSeq . ";size=" . $row->{count} . ";\n" . $row->{sequence} . "\n";
      $nSeq++;
    }
    close($fhFna) || die "Error writing to $tmpFna\n";
    my $usearchCmd = "$usearch -unoise3 $tmpFna -minsize $minN -ampout $tmpU";
    $usearchCmd .= " -quiet" unless defined $debug;
    system($usearchCmd) == 0 || die "$usearchCmd -- failed: $!";
    open(my $fhU, "<", $tmpU) || die "Cannot read $tmpU\n";
    while (my $line = <$fhU>) {
      next unless $line =~ m/^>/;
      chomp $line;
      die "Cannot parse sequence number and amptype from -ampout file $tmpU line $line"
        unless $line =~ m/^>SEQ(\d+);size=\d+;amptype=([a-z]+)/i;
      my $num = $1;
      my $amptype = $2;
      die "Unknown sequence number $num" unless exists $numToRow{$num};
      next unless $amptype eq "otu"; # ignore chimeras
      my $row = $numToRow{$num};
      if (!exists $seqToEsv{$row->{sequence}}) {
        my $newName;
        if ($namePre ne "") {
          $newName = $namePre . (1 + scalar(keys %seqToEsv));
        } else {
          $newName = md5_hex($row->{sequence});
        }
        die "zotu name $newName is already in use!"
          if exists $esvToSeq{$newName};
        $seqToEsv{$row->{sequence}} = $newName;
        $esvToSeq{$newName} = $row->{sequence};
      }
      my $esv = $seqToEsv{ $row->{sequence} };
      die unless defined $esv;
      die if exists $counts{ $row->{primer_name} }{ $index }{ $esv };
      $counts{ $row->{primer_name} }{ $index }{ $esv } = $row->{count};
    }
    close($fhU) || die "Error reading $tmpU";
  }
}

unlink($tmpU);
unlink($tmpFna);

print STDERR "Found a total of " . scalar(keys %esvToSeq) . " ESVs across " . scalar(@infiles) . " input files\n";
open(my $fhFna, ">", "$outPre.fna")
  || die "Cannot write to $outPre.fna\n";
foreach my $esv (sort compareEsvNames keys %esvToSeq) {
  die $esv unless $esvToSeq{$esv};
  print $fhFna ">" . $esv . "\n" . $esvToSeq{$esv} . "\n";
}
close($fhFna) || die "Error writing to $outPre.fna";
print STDERR "Wrote $outPre.fna\n";

open(my $fhTab, ">", "$outPre.tsv")
  || die "Cannot write to $outPre.tsv";
print $fhTab join("\t", "primer_name", "index", "Zotu", "count", "tot")."\n";
foreach my $primer (sort keys %counts) {
  my $hashP = $counts{$primer};
  foreach my $index (sort keys %$hashP) {
    my $hashI = $counts{$primer}{$index};
    my $tot = sum(values %$hashI);
    $totReadsOut += $tot;
    foreach my $esv (sort compareEsvNames keys %$hashI) {
      print $fhTab join("\t", $primer, $index, $esv, $hashI->{$esv}, $tot)."\n";
    }
  }
}
close($fhTab) || die "Error writing to $outPre.tab";
print STDERR "Wrote $outPre.tsv\n";
print STDERR sprintf("Total kept count: %.2fM of %.2fM input (%.1f%%)\n",
                     $totReadsOut/1e6, $totReadsIn/1e6, 100*($totReadsOut+1)/($totReadsIn+1));

sub compareEsvNames($$) {
  my ($a,$b) = @_;
  my $aIsPre = substr($a, 0, length($namePre)) eq $namePre;
  my $bIsPre = substr($b, 0, length($namePre)) eq $namePre;
  if ($aIsPre && $bIsPre) {
    my $a2 = substr($a, length($namePre));
    my $b2 = substr($b, length($namePre));
    return $a2 <=> $b2 if $a2 =~ m/^\d+$/ && $b2 =~ m/^\d+$/;
  }
  #else
  return $a cmp $b;
}
