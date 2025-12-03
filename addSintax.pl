#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib $RealBin;
use dmUtils qw{ReadTable ReadColumnNames};

my $usearch = "$RealBin/usearch";
my $db = "$RealBin/rdp_16s_v18.udb";

my $usage = <<END
Usage: addSintax.pl -in clean [ -out clean.sintax.tsv ]

  Given clean.fna and clean.tsv from parseInline.pl, use sintax (from
  usearch) to add taxonomic assignments to the table, with the fields
  "taxLong" (with confidence values) and "domain", "phylum", "class",
  "order", "family", and "genus" (with only the high-confidence
  assignments included). Only sintax hits on the - strand are considered.

Additional options:
  -usearch $usearch
     the usearch executable
  -db $db
     the sintax database
END
  ;

my ($in, $outFile);
die $usage
  unless GetOptions('in=s' => \$in,
                    'out=s' => \$outFile,
                    'usearch=s' => \$usearch,
                    'db=s' => \$db)
  && @ARGV == 0;
die $usage unless defined $in;
my $inTable = "$in.tsv";
my $inFasta = "$in.fna";
foreach my $file ($inTable,$inFasta,$db) {
  die "No such file: $file\n" unless -e $file;
}
die "No such executable: $usearch\n"
  unless -x $usearch;

$outFile = "$in.tax.tsv" unless defined $outFile;

my @in = ReadTable($inTable, ["Zotu"]);

my $tmpFile = "/tmp/addSintax.$$.tsv";
my @cmd = ($usearch, "-sintax", $inFasta, "-db", $db, "-strand", "both",
       "-tabbedout", $tmpFile,
       "-quiet");
system(@cmd) == 0
  || die "Error running @cmd -- $!";

my %tax = (); # zotu => [ taxLong, tax ]
my %skip = (); # zotu => 1 if has no hit or on the wrong strand, and is in the input table
open(my $fh, "<", $tmpFile) || die "Cannot read $tmpFile\n";
while (my $line = <$fh>) {
  chomp $line;
  my ($zotu, $taxLong, $strand, $tax) = split /\t/, $line;
  die "Invalid line from usearch: $line" unless defined $tax;
  if ($strand eq "-") {
    die "Duplicate zotu $zotu in $tmpFile" if exists $tax{$zotu};
    $tax{$zotu} = [ $taxLong, $tax ];
  }
}
close($fh) || die "Error reading $tmpFile\n";
unlink($tmpFile);

my @inFields = ReadColumnNames($inTable);
my @taxPartNames = qw{domain phylum class order family genus};
my @taxPartKeys = map substr($_, 0, 1), @taxPartNames;
open (my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
print $fhOut join("\t", @inFields, @taxPartNames, "taxLong")."\n";
foreach my $row (@in) {
  my $zotu = $row->{"Zotu"};
  my ($taxLong, $tax);
  if (exists $tax{$zotu}) {
    ($taxLong, $tax) = @{ $tax{$zotu} };
  } else {
    ($taxLong, $tax) = ("", "");
    $skip{$zotu} = 1;
  }
  my @taxParts = split /,/, $tax;
  my %taxParts = ();
  foreach my $part (@taxParts) {
    $part =~ m/^([a-z]):(.*)$/ || die "Cannot parse tax part $part from $tax";
    my ($key, $value) = ($1, $2);
    $value =~ s/^"//;
    $value =~ s/"$//;
    $taxParts{$key} = $value;
  }
  my @taxShow = map { $taxParts{$_} || "" } @taxPartKeys;
  my @inValues = map $row->{$_}, @inFields;
  print $fhOut join("\t", @inValues, @taxShow, $taxLong)."\n";
}
close($fhOut) || die "Error writing to $outFile\n";
print STDERR "Wrote $outFile -- no tax assignment for " . scalar(keys %skip) . " Zotus\n";
