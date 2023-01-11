#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
sub Run($$);
sub FileToSample($);
sub IsFileEmpty($); # checks for empty gzipped files as well (which do not have size 0)

my $pear = "$RealBin/pear";
my $usearch = "$RealBin/usearch";
my $parser = "$RealBin/parseInline.pl";
my $submitter = "$RealBin/submitter.pl";
my $maxE = 1;

my $usage = <<END
Usage: runInline.pl -model 806R -dir indir
  Given a directory with paired-end amplicon reads, builds tables of
  the sequences seen in each sample.  Uses pear to merge reads,
  usearch to filter out noisy reads, and parseInline.pl to tabulate
  the sequences seen in each sample.

Optional arguments:
-test -- report what commands to run (to stdout), but do no work
-maxE -- maximum errors per read (default $maxE)
-modelFile -- specify the model by file instead of by name

You can alter the executables used with:
-pear $pear
-usearch $usearch
-parser $parser

You can also set -endSeq and -endRange -- see parseInline.pl
END
;

my ($dir, $modelName, $modelFile, $endSeq, $endRange, $test);
die $usage
  unless GetOptions('dir=s' => \$dir,
                    'model=s' =>\$modelName,
                    'modelFile=s' => \$modelFile,
                    'test' => \$test,
                    'maxE=f' => \$maxE,
                    'pear=s' => \$pear,
                    'usearch=s' => \$usearch,
                    'parser=s' => \$parser,
                    'endSeq=s' => \$endSeq,
                    'endRange=s' => \$endRange)
  && defined $dir;
die "Invalid max errors, must be positive\n" if $maxE <= 0;
die "Not a directory: $dir\n" unless -d $dir;
die "Must specify -model or -modelFile (but not both)\n"
  unless (defined $modelName) xor (defined $modelFile);

unless (defined $modelFile) {
  $modelFile = "inline_${modelName}.tsv";
  if (! -e $modelFile) {
    my $modelFileOrig = $modelFile;
    $modelFile = "$RealBin/$modelFile";
    die "No such file: $modelFileOrig or $modelFile\n"
      unless -e $modelFile;
  }
}
die "No such file: $modelFile\n" unless -e $modelFile;

foreach my $x ($pear,$usearch,$parser,$submitter) {
  die "No such executable: $x\n" unless -x $x;
}

my @fastq = glob("$dir/*.fastq.gz");
if (@fastq == 0) {
  @fastq = grep !m#[.]pear[.][^/]+$#, glob("$dir/*.fastq");
}
die "No fastq inputs found in $dir\n" if @fastq == 0;
@fastq = sort @fastq;
@fastq = grep !m#/Undetermined[^/]+$#i, @fastq;

my @r1 = grep m/R1_\d+[.]fastq[.gz]+$/, @fastq;
die "No R1 input reads found in $dir\n" if @r1 == 0;
my @r2 = grep m/R2_\d+[.]fastq[.gz]+$/, @fastq;

die "Found different numbers of fastq files for R1 and R2\n"
  unless scalar(@r1) == scalar(@r2);

print STDERR "Found " . scalar(@r1) . " input files\n";

my @samples = (); # a sample name from each read
for (my $i = 0; $i < scalar(@r1); $i++) {
  my $r1 = $r1[$i];
  my $r2 = $r2[$i];
  my $s1 = FileToSample($r1);
  my $s2 = FileToSample($r2);
  die "Mismatched read files, $r1 and $r2\n" unless $s1 eq $s2;
    push @samples, $s1;
}

my %samples = map { $_ => 1 } @samples;
die "Not all sample identifiers are unique. Multiple read files for one sample are not yet supported.\n"
  unless scalar(keys %samples) == scalar(@samples);
print STDERR "Sample names: " . join(" ", @samples) . "\n";

my @pearCmds = ();
my @usearchCmds = ();
my @parseCmds = ();
for(my $i = 0; $i < scalar(@r1); $i++) {
  my $r1 = $r1[$i];
  my $r2 = $r2[$i];
  my $sample = $samples[$i];
  if (IsFileEmpty($r1)) {
    die "File $r1 is empty but $r2 is not\n" unless IsFileEmpty($r2);
    print STDERR "Skipping sample $sample with no reads\n";
  } else {
    push @pearCmds, "$pear -f $r1 -r $r2 -o $dir/$sample.pear >& $dir/$sample.pear.log";
    push @usearchCmds, "$usearch -fastq_filter $dir/$sample.pear.assembled.fastq -fastq_maxee $maxE -fastaout $dir/$sample.filtered";
    push @parseCmds, "($parser -modelFile $modelFile < $dir/$sample.filtered > $dir/$sample.parse.tab) >& $dir/$sample.parse.log";
  }
}

if (defined $test) {
  print "# pear commands\n";
  print join("\n", @pearCmds) . "\n";
  print " usearch commands\n";
  print join("\n", @usearchCmds) . "\n";
  print "# parse commands\n";
  print join("\n", @parseCmds)."\n";
} else {
  Run(\@pearCmds, "$dir/pear.cmds");
  Run(\@usearchCmds, "$dir/usearch.cmds");
  Run(\@parseCmds, "$dir/parse.cmds");
  print STDERR "Finished demultiplexing in $dir -- see $dir/*.parse.tab\n";
  print STDERR "To remove noisy or chimeric reads, use cleanInline.pl, i.e.\n";
  print STDERR "cleanInline.pl -in $dir/*.parse.tab -out $dir/clean\n";
}

sub Run($$) {
  my ($cmds, $cmdsFile) = @_;
  open(my $fh, ">", $cmdsFile) || die "Cannot write to $cmdsFile\n";
  print $fh join("\n", @$cmds)."\n";
  close($fh) || die "Error writing to $cmdsFile\n";
  system($submitter, $cmdsFile) == 0
    || die "Error running commands in $cmdsFile\n";
}

# Except sample names of the form IT001 or B09
sub FileToSample($) {
  my ($fileName) = @_;
  my $name = $fileName;
  $name =~ s!^.*/!!;
  $name =~ s/_R[12]_\d+//;

  die "Cannot extract sample number from $name"
    unless $name =~ m/^(IT\d+)/ || $name =~ m/[_-](IT\d+)/
      || $name =~ m/^(S\d+)/ || $name =~ m/_(S\d+)/
      || $name =~ m/^([A-Z][0-9]+)[_.]/;
  die "Cannot parse sample name from file $fileName name $name" unless $1;
  return $1;
}

sub IsFileEmpty($) {
  my ($file) = @_;
  my $fh;
  if ($file =~ m/[.]gz$/) {
    open($fh, "zcat $file |") || die "Cannot gunzip $file\n";
  } else {
    open($fh, "<", $file) || die "Cannot read $file\n";
  }
  my $line = <$fh>;
  close($fh);
  return $line ? 0 : 1;
}
