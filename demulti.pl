#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib $RealBin;
use dmUtils qw{ReadTable reverseComplement};

my $endSeq = "TTACCGCGGCKGCTGRCAC";
my $endRange = "1:4";

my $usage = <<END
Usage: demulti.pl -model 926R -reads read1file read2file  -expect 1,2,3,8 -out prefix

The model must correspond to a file named inline_model.tsv (in either
the current directory or the script directory) that includes the
fields primer_name, Ns, inline_index, and begin. Usually it is either
806R or 926R.

There must be two input fastq files, one for each read. If their names
end with .gz they are assumed to be compressed.

-expect indicates which primers are expected. Only output files for
 these numbers will be created. I.e. -expect 1 refers to the primer
 whose name ends with 1.

The output files for primer 1 will be named out_p1_R1.fastq and out_p1_R2.fastq

Optional arguments:
-end $endSeq -- the end sequence  (whose reverse complement
  should be at the beginning of read 2)
-end $endRange -- number of Ns after the end sequence
-noTrimEnd -- skip trimming of the end sequence
-debug --verbose output to standard error for every read
END
;

my @readFiles;
my ($modelFile, $modelName, $outPrefix, $noTrimEnd, $debug, $expectSpec);
die $usage
  unless @ARGV > 0
  && GetOptions('model=s' => \$modelName,
                'modelFile=s' => \$modelFile,
                'end=s' => \$endSeq,
                'endRage=s' => \$endRange,
                'debug' => \$debug,
                'expect=s' => \$expectSpec,
                'noTrimEnd' => \$noTrimEnd,
                'reads=s{2,2}' => \@readFiles,
                'out=s' => \$outPrefix)
  && @ARGV == 0;
die "Must specify -model or -modelFile (but not both)\n"
  unless (defined $modelName) xor (defined $modelFile);
die "Invalid end range $endRange: should be something like 1:4\n" unless $endRange =~ m/^\d+:\d+$/;
my ($minEnd,$maxEnd) = split /:/, $endRange;
die "Invalid end range $endRange" unless $maxEnd >= $minEnd;
die "Invalid end sequence $endSeq\n" unless $endSeq =~ m/^[A-Z]+$/;
my $endPattern = $endSeq;
$endPattern =~ s/[^ACGT]/./g; # treat all ambiguity characters as any character

my $endPatternRev = reverseComplement($endPattern);

die "Invalid expect specifier $expectSpec\n" unless $expectSpec =~ m/^[0-9,]+$/;
my @expectSpec = split /,/, $expectSpec;
my %expectSpec = map { $_ => 1 } @expectSpec;

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
my %inlineSeq = (); # number => object
foreach my $inline (@inline) {
  my $num = $inline->{primer_name};
  $num =~ s/^.*[a-zA-Z_](\d+)$/$1/ || die "Invalid primer name $num has no number";
  die "duplicate inline number $num" if exists $inlineSeq{$num};
  $inline->{num} = $num;
  $inlineSeq{$num} = $inline;
}
my %primers = (); # primer name to row
my $inlineLen;
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

foreach my $expect (@expectSpec) {
  die "Expect number $expect does not appear in the input table $modelFile\n"
    unless exists $inlineSeq{$expect};
}

foreach my $readFile (@readFiles) {
  die "No such file: $readFile\n"
    unless -e $readFile;
}



print STDERR "Parsing reads for primers $expectSpec from\n  $readFiles[0]\n  $readFiles[1]\ninto ${outPrefix}_p*_R*.fastq\n";
print STDERR "Trimming end $endSeq (rc $endPatternRev) with N${minEnd}:${maxEnd}\n" unless defined $noTrimEnd;

my ($fhFwd, $fhRev);
if ($readFiles[0] =~ m/[.]gz$/) {
  open($fhFwd, "-|", "zcat", $readFiles[0]) || die "Cannot zcat $readFiles[0]\n";
} else {
  open($fhFwd, "<", $readFiles[0]) || die "Cannot read $readFiles[0]\n";
}
if ($readFiles[1] =~ m/[.]gz$/) {
  open($fhRev, "-|", "zcat", $readFiles[1]) || die "Cannot zcat $readFiles[1]\n";
} else {
  open($fhRev, "<", $readFiles[1]) || die "Cannot read $readFiles[1]\n";
}

my %fhOut = (); # sample number to [fwdOut, revOut]
my %fileOut = (); # sample number to [fwdFile,revFile]
foreach my $exp (@expectSpec) {
  foreach my $name ("R1","R2") {
    my $file = "${outPrefix}_p${exp}_${name}.fastq";
    open(my $fh, ">", $file) || die "Cannot write to $file\n";
    push @{ $fileOut{$exp} }, $file;
    push @{ $fhOut{$exp} }, $fh;
  }
}

my $nRead = 0;
my $nKept = 0;
my $nNoEnd = 0;
while(my $headerFwd = <$fhFwd>) {
  $nRead++;
  chomp $headerFwd;
  die "Header does not start with @ -- $headerFwd" unless $headerFwd =~ m/^@/;
  my $headerRev = <$fhRev> || die "No corresponding 2nd read for $headerFwd";
  chomp $headerRev;
  my $short = $headerFwd; $short =~ s/ .*//;
  my $short2 = $headerRev; $short2 =~ s/ .*//;
  die "Mismatching headers: $short vs. $short2" unless $short eq $short2;
  my $seqFwd = <$fhFwd> || die "Missing sequence in $readFiles[0]";
  my $seqRev = <$fhRev> || die "Missing sequence in $readFiles[1]";
  my $skipFwd = <$fhFwd> || die "Missing line in $readFiles[0]";
  my $skipRev = <$fhRev> || die "Missing line in $readFiles[1]";
  my $qFwd = <$fhFwd> || die "Missing quality in $readFiles[0]";
  my $qRev = <$fhRev> || die "Missing quality in $readFiles[0]";
  chomp $seqFwd;
  chomp $seqRev;
  chomp $skipFwd;
  chomp $skipRev;
  chomp $qFwd;
  chomp $qRev;
  die "Wrong length for quality in $readFiles[0]" unless length($seqFwd) == length($qFwd);
  die "Wrong length for quality in $readFiles[1]" unless length($seqRev) == length($qRev);

  # Classify the sequence using the beginning of $seqFwd
  my ($thisPrimer, $cutTo);
  foreach my $at (0..(scalar(@inlineSeq)-1)) {
    my $beginHash = $inlineSeq[$at];
    next unless defined $beginHash;
    my $beginThis = substr($seqFwd, $at, $inlineLen);
    if (exists $beginHash->{$beginThis}) {
      my $primerName = $beginHash->{$beginThis};
      my $beginLen = length($primers{$primerName}{begin});
      my $begin2 = substr($seqFwd, $at + $inlineLen, $beginLen);
      if ($begin2 =~ $primers{ $primerName }{beginPattern}) {
        print STDERR "Demultiplex $headerFwd into $primerName -- starts with "
          . substr($seqFwd, 0, $at) . " $beginThis $begin2\n"
          if $debug;
        $thisPrimer = $primerName;
        $cutTo = $at + $inlineLen + $beginLen;
        last;
        next;
      }
    }
  }

  if (defined $thisPrimer) {
    my $num = $thisPrimer;
    $num =~ s/^.*[a-zA-Z_](\d+)$/$1/ || die "Invalid primer name $num has no number";
    my $revCutTo;
    if (exists $fhOut{$num}) {
      # Truncate the end
      unless (defined $noTrimEnd) {
        my $revCheck = substr($seqRev, 0, length($endSeq) + $maxEnd);
        if ($revCheck =~ $endPatternRev) {
          $revCutTo = $-[0] + length($endSeq);
        } else {
          $nNoEnd++;
          print STDERR "No end sequence for $headerRev\n" if defined $debug;
        }
      }
      if (defined $revCutTo || defined $noTrimEnd) {
        $nKept++;
        my ($fhFwd,$fhRev) = @{ $fhOut{$num} };
        $seqFwd = substr($seqFwd, $cutTo);
        $qFwd = substr($qFwd, $cutTo);
        if (defined $revCutTo) {
          $seqRev = substr($seqRev, $revCutTo);
          $qRev = substr($qRev, $revCutTo);
        }
        print $fhFwd $headerFwd, "\n", $seqFwd, "\n", $skipFwd, "\n", $qFwd, "\n";
        print $fhRev $headerRev, "\n", $seqRev, "\n", $skipRev, "\n", $qRev, "\n";
      }
    } elsif (defined $debug) { # primer is unexpected
      print STDERR "Skipped unexpected demultiplex for $thisPrimer\n";
    }
  } elsif (defined $debug) { # not demultiplexed
    print STDERR "Failed to demultiplex $headerFwd beginning with " . substr($seqFwd, 0, 20) . "\n";
  }
}
close($fhFwd) || die "Error reading $readFiles[0]";
close($fhRev) || die "Error reading $readFiles[1]";
foreach my $exp (@expectSpec) {
  foreach my $fh (@{ $fhOut{$exp} }) {
    close($fh) || die "Error writing output files for $exp";
  }
}
print STDERR "Read $nRead reads, kept $nKept\n";
print STDERR "Multiplexed but missing end for $nNoEnd\n" unless defined $noTrimEnd;
