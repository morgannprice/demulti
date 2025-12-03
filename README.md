# Scripts for dual demultiplexing of 16S amplicon reads

These scripts are intended to demultiplex paired-end MiSeq reads with
independent indexes on both sides. Because of index hopping, using
these primers with more recent Illumina platforms is not recommended.

# System Requirements

These scripts should work on any Linux system, and would probably work
on other Unix or MacOS as well. All of the code is written in perl. (I
use perl v5.16.3.)

If you just want to demultiplex the paired-end reads into new fastq
files (demulti.pl), then there are no other requirements. If you want
to identify exact sequence variants, then you'll also need to install
[PEAR](https://cme.h-its.org/exelixis/web/software/pear/) and
[usearch](http://www.drive5.com/usearch/). You should put these
executables, or symbolic links to these executables, in the demulti/
directory.

(I use PEAR 0.9.11 and usearch v10.0.240_i86linux32.)

# Input format

The paired-end reads must be in fastq format or gzipped fastq
format. For runInline.pl, the reads must have names ending in
R1_001.fastq.gz and R2_001.fastq.gz. Other numeric suffixes instead of
001 are allowed. Omit .gz if the files are not compressed.

# Primers for V4 or V4V5 

The primers corresponding to the 5' end of the reads are described in
inline_806R.tsv (for V4 only) or inline_926R.tsv (for V4V5). They
include the "P5" adapter for Illumina sequencing, 1-4 Ns, 8
nucleotides that vary by the index (inline_index), and an initial
sequence that matches to 16S sequences (primer_seq). The portion of
the primer that matches the 16S is either GGACTACNVGGGTWTCTAAT (806R)
or CCGYCAATTYMTTTRAGTTT (926R).

The primers corresponding to the 3' end of the reads are described in
515F_primers.tsv. They include the P7 adapter, the 6-nucleotide
Illumina index (IT001:IT096), a sequence for priming the index read,
1-4 Ns, and an initial sequence that matches to 16S at 515F
(GTGYCAGCMGCCGCGGTAA). These are similar to the updated primers from
the
[Earth Microbiome Project](https://earthmicrobiome.org/protocols-and-standards/16s/).

These primer pairs were designed by Adam Deutschbauer and Hans Carlson.

# Scripts

There are two independent ways to analyze the reads:

* Demultiplex to paired-end fastq with demulti.pl
* Demultiplex and identify exact sequence variants with runInline.pl and cleanInline.pl

If you identify exact sequence variants, there's also a script to add
taxonomic classification from sintax (addSintax.pl).

# Demultiplexing to fastq with demulti.pl

Run demulti.pl for each input file (i.e., IT001-IT096) and specify which inline primers are expected. This will not filter the reads in any way -- every sequence that has the expecting beginning and ending sequences (including exact matches to the primers) is extracted. (This could produce a large number of fastq files.) Only the portion of the reads inside the primers is maintained.

To analyze all 96 input files, you can use a bash for loop like this to analyze files with names like IT001_S1_L001_R1_001.fastq.gz:

```
(for i in `seq 1 96`; do
nice demulti.pl -model 926R -reads IT*${i}_S${i}_L001_R{1,2}_001.fastq.gz -expect 1,3 -out demulti_S${i};
done) >& demulti.log
```

# Identifying exact sequence variants with runInline.pl and cleanInline.pl

These scripts will remove sequences that are likely to have errors,
sequences that are rare, and sequences that are likely to be
chimeras. The defaults are: filter out sequences with more than 1
expected error and that occur less than 4 times in a sample. Likely chimeras are removed using
[unoise3](https://www.drive5.com/usearch/manual/unoise_algo.html).

First use runInline.pl, then use cleanInline.pl. runInline.pl expects
that specified directory (-dir) includes fastq or fastq.gz files and
that the file names include words like IT001 or A01. runInline.pl
takes a while, so you might want to run it in the background. For
cleanInline.pl, tell it which inline indexes you expect to be present. Example:

```
nice runInline.pl -dir . -model 926R >& inline.log
cleanInline.pl -in ./*parse.tab -out mytable -primers 23,24
```

This will produce a fasta file named mytable.fna with all the exact
sequence variants, and a tab-delimited table named mytable.tsv with
how often each sequence was observed in each sample.

# Adding taxonomic information with addSintax.pl

Once you have the output of cleanInline.pl, you can use addSintax.pl to add taxonomic information from usearch's
[sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html).

For addSintax.pl to work, you need to make the sintax version of the
Ribosomal Database Project, which you can do with these commands:

```
wget https://www.drive5.com/sintax/rdp_16s_v18.fa.gz
gunzip rdp_16s_v18.fa.gz
usearch -makeudb_sintax rdp_16s_v18.fa -output rdp_16s_v18.udb
```

Once the sintax database is built, you can add taxonomic information with

```
addSintax.pl -in mytable
```

# Contact

If you find any bugs in these scripts, please let me know.
You can reach me at funwithwords26@gmail.com. Thanks, Morgan.
