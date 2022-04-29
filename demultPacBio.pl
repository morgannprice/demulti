#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib $RealBin;
use dmUtils qw{ReadFastaEntry reverseComplement};

sub listVariants($); # return all 1-nt variants for a sequence
sub matchBarcode($$$);

my $bcLen = 16;
my $slop = 2;
my $bcFile = "$RealBin/barcodesPacBio.fna";
my $lenExpect = 1453;
my $lenRange = 300;

my $leftPattern = "AGnGTTnGATnnTGGCTCAG";
my $rightPattern = "AAGTCGTAACAAGGTAnC";

my $usage = <<END
demultPacBio.pl < filtered.fna > parsed.tsv
  Given a fasta file of CCS reads, find the 16S flanking sequences and barcodes,
    and makes a table of inserts (the 16S region within the primers) and barcodes.
  Also tries to extract expected errors (ee) from the header.
  Ignores reads that do not have the expected 16S flanking sequences or insert size.
  Reads are checked in both orientations, and reads may have any amount of
    other sequences on the sides.
  Barcodes are allowed to have up to 1 error or to be truncated by up to 2 nt.
Optional arugments:
-bc barcodeFile -- default $bcFile
-bcLen $bcLen -- barcode size
-slop $slop -- how far from expected location to check for barcodes.
-len $lenExpect -- expected length of insert
-range $lenRange -- how much length is allowed to vary by
-left $leftPattern -- expected pattern at left
-right $rightPattern -- expected pattern at right
END
;

die $usage
  unless GetOptions('bcLen=i' => \$bcLen,
                    'slop=i' => \$slop,
                    'bcFile=s' => \$bcFile,
                    'len=i' => \$lenExpect,
                    'range=i' => \$lenRange,
                    'left=s' => \$leftPattern,
                    'right=s' => \$rightPattern)
  && @ARGV == 0;
$leftPattern =~ s/[^ACGT]/./g;
$rightPattern =~ s/[^ACGT]/./g;

die "No file: $bcFile\n" unless -e $bcFile;

my %bc = (); # seq => name
open(my $fh, "<", $bcFile) || die "Cannot read $bcFile\n";
my $state = {};
while (my ($header, $seq) = ReadFastaEntry($fh, $state)) {
  die "barcode $header is not the expected length\n"
    unless length($seq) == $bcLen;
  my $rc = reverseComplement($seq);
  $bc{$seq} = $header;
  $bc{$rc} = $header;
  # allow up to  2 nt at either end to be missing
  $bc{substr($seq, 0, $bcLen-2)} = $header;
  $bc{substr($seq, 2)} = $header;
  $bc{substr($rc, 0, $bcLen-2)} = $header;
  $bc{substr($rc, 2)} = $header;
  # allow all single-nucleotide errors
  foreach my $v (listVariants($seq)) {
    $bc{$v} = $header;
  }
  foreach my $v (listVariants($rc)) {
    $bc{$v} = $header;
  }
}
close($fh) || die "Error reading $bcFile\n";

print STDERR "Looking for $leftPattern and $rightPattern\n";

print join("\t", qw{read ee strand barcode1 barcode2 readLen insertAt len seq})."\n";
$state = {};
my $nSkip = 0;
my $nRead = 0;
while (my ($header, $seq) = ReadFastaEntry(\*STDIN, $state)) {
  $nRead++;
  my $ee = "";
  if ($header =~ m/^(.*);ee=(.*)$/) {
    ($header, $ee) = ($1,$2);
  }
  my $strand = "+";
  unless ($seq =~ m/$leftPattern/ && $seq =~ m/$rightPattern/) {
    $seq = reverseComplement($seq);
    $strand = "-";
  }
  $seq =~ m/$leftPattern/g;
  my $leftAt = pos $seq;
  $seq =~ m/$rightPattern/g;
  my $rightAt = pos $seq;

  unless (defined $leftAt && defined $rightAt && $rightAt > $leftAt + length($rightPattern)) {
    $nSkip++;
    next;
  }

  my $insert = substr($seq, $leftAt, $rightAt - $leftAt - length($rightPattern));
  unless (length($insert) >= $lenExpect - $lenRange
          && length($insert) <= $lenExpect + $lenRange) {
    $nSkip++;
    next;
  }
  my $left = substr($seq, 0, $leftAt - length($leftPattern));
  my $right = substr($seq, $rightAt);

  my $leftPos = length($left) - $bcLen - $slop;
  $leftPos = 0 if $leftPos < 0;
  my $barcodeL = matchBarcode(substr($left, $leftPos), \%bc, $slop);
  my $barcodeR = matchBarcode(substr($right, 0, $bcLen + $slop), \%bc, $slop);
  print join("\t", $header, $ee, $strand, $barcodeR, $barcodeL,
             length($seq), $leftAt, length($insert), $insert)."\n";
}
print STDERR "Skipped $nSkip of $nRead\n";


sub matchBarcode($$$) {
  my ($seq, $bcHash, $slop) = @_;
  return "" unless length($seq) >= $bcLen - 2;
  for (my $at = 0; $at <= $slop; $at++) {
    my $bcSeq = substr($seq, $at, $bcLen);
    my @bcTest = ($bcSeq, substr($bcSeq, 2), substr($bcSeq, 0, $bcLen-2));
    foreach my $bcTest (@bcTest) {
      return $bcHash->{$bcTest} if exists $bcHash->{$bcTest};
    }
  }
  return "";
}

sub listVariants($) {
    my ($baseseq) = @_;
    my @out = ();
    $baseseq = uc($baseseq);
    my $len = length($baseseq);
    foreach my $i (0..($len-1)) {
        my $pre = substr($baseseq,0,$i);
        my $char = substr($baseseq,$i,1);
        my $post = substr($baseseq,$i+1);
        next unless $char eq "A" || $char eq "C" || $char eq "G" || $char eq "T";
        foreach my $newchar (qw{A C G T}) {
            push @out, $pre . $newchar . $post unless $newchar eq $char;
        }
    }
    return(@out);
}
