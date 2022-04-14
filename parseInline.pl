#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib $RealBin;
use dmUtils qw{ReadTable ReadFastaEntry};

my $endSeq = "TTACCGCGGCKGCTGRCAC";
my $endRange = "1:4";
my $debug;

my $usage = <<END
Usage: parseInline.pl -model 806R < filtered_fasta > output.tsv

The model must correspond to a file named inline_model.tsv (in either
the current directory or the script directory) that includes the fields
primer_name, Ns, inline_index, and begin

Optional arguments:
-end $endSeq
-endRange $endRange -- number of Ns after the end sequence
-modelFile filename -- this can be used instead of -model
-debug -- verbose output to standard error for every read
END
;

my ($modelFile, $modelName);
die $usage
  unless @ARGV > 0
  && GetOptions('model=s' => \$modelName,
                'modelFile=s' => \$modelFile,
                'end=s' => \$endSeq,
                'endRange=s' => \$endRange,
                'debug' => \$debug)
  && @ARGV == 0;
die "Must specify -model or -modelFile (but not both)\n"
  unless (defined $modelName) xor (defined $modelFile);
die "Invalid end range $endRange: should be something like 1:4\n" unless $endRange =~ m/^\d+:\d+$/;
my ($minEnd,$maxEnd) = split /:/, $endRange;
die "Invalid end range $endRange" unless $maxEnd >= $minEnd;
die "Invalid end sequence $endSeq\n" unless $endSeq =~ m/^[A-Z]+$/;
my $endPattern = $endSeq;
$endPattern =~ s/[^ACGT]/./g; # treat all ambiguity characters as any character

unless (defined $modelFile) {
  $modelFile = "inline_${modelName}.tsv";
  if (! -e $modelFile) {
    my $modelFileOrig = $modelFile;
    $modelFile = "$RealBin/$modelFile";
    die "No such file: $modelFileOrig or $modelFile\n"
      unless -e $modelFile;
  }
}

my @inline = ReadTable($modelFile, ["primer_name", "Ns", "inline_index", "begin"]);
my @inlineSeq = (); # position (0-based) => expected sequence => primer name
my $inlineLen;
my %primers = (); # primer name to row
foreach my $row (@inline) {
  die "Invalid Ns specifier $row->{Ns}" unless $row->{Ns} =~ m/^N+$/i;
  my $at = length($row->{Ns});
  die "Invalid index $row->{inline_index}" unless $row->{inline_index} =~ m/^[ACGT]+$/;
  die "inline_index must all be the same length"
    if defined $inlineLen && length($row->{inline_index}) != $inlineLen;
  $inlineLen = length($row->{inline_index});
  die "Invalid primer_name $row->{primer_name}" if $row->{primer_name} eq "";
  die "Duplicate primer_name $row->{primer_name}" if exists $primers{ $row->{primer_name} };
  die "Invalid begin $row->{begin}" unless $row->{begin} =~ m/^[A-Z]+$/;
  $row->{beginPattern} = $row->{begin};
  $row->{beginPattern} =~ s/[^ACGT]/./g;
  $primers{ $row->{primer_name} } = $row;
  $inlineSeq[$at]{ $row->{inline_index} } = $row->{primer_name};
}

my %counts = (); # primer name => sequence => count

my $nReads = 0;
my $nSkip = 0;
my $nEndFail = 0;
my $state = {};
while (my ($header, $seq) = ReadFastaEntry(\*STDIN, $state)) {
  if ($seq eq "") { # occurs sometimes in testing shortened files
    print STDERR "Warning! Ignoring empty sequence for $header\n";
    next;
  }
  my $use = 0;
  # Classify the sequence
  foreach my $at (0..(scalar(@inlineSeq)-1)) {
    my $beginHash = $inlineSeq[$at];
    next unless defined $beginHash;
    my $beginThis = substr($seq, $at, $inlineLen);
    if (exists $beginHash->{$beginThis}) {
      my $primerName = $beginHash->{$beginThis};
      my $beginLen = length($primers{$primerName}{begin});
      my $begin2 = substr($seq, $at + $inlineLen, $beginLen);
      if ($begin2 =~ $primers{ $primerName }{beginPattern}) {
        print STDERR "Demultiplex $header into $primerName -- starts with "
          . substr($seq, 0, $at) . " $beginThis $begin2\n"
          if $debug;
        my $remain = substr($seq, $at + $inlineLen + $beginLen);
        my $endLen = $maxEnd + length($endSeq);
        my $endThis = substr($remain, length($remain) - $maxEnd - length($endSeq));
        if ($endThis =~ $endPattern) {
          my $remain2 = substr($remain, 0, length($remain) - $maxEnd - length($endSeq) + $-[0]);
          $counts{$primerName}{$remain2}++;
          $use = 1;
          print STDERR "Using $header in $primerName -- ends with $endThis; keeping "
            . substr($remain2, 0, 10) . "..." . substr($remain2, length($remain2)-10, 10) . "\n"
            if $debug;
          last;
        } else {
          print STDERR "End fail for $header -- no match for $endPattern in $endThis\n"
            if $debug;
          $nEndFail++;
          last; # confirmed this was the correct primer even though not used
        }
      }
    }
  }
  $nReads++;
  $nSkip++ unless $use;
  print STDERR "Skipping $header starts with " . substr($seq, 0, $inlineLen + scalar(@inlineSeq) + 10) . "\n"
    if $debug && ! $use;
}
print STDERR "Ignoring $nSkip of $nReads reads\n$nEndFail of the $nSkip demultiplexed but did not match the expected end\n";

print join("\t", "primer_name", "count", "sequence")."\n";
foreach my $row (@inline) {
  my $primerName = $row->{primer_name};
  next unless exists $counts{$primerName};
  my $countHash = $counts{$primerName};
  my @seqs = sort keys %$countHash;
  @seqs = sort { $countHash->{$b} <=> $countHash->{$a} } @seqs;
  foreach my $seq (@seqs) {
    print join("\t", $primerName, $countHash->{$seq}, $seq)."\n";
  }
}
