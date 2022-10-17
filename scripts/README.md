# oNTORRA Scripts

## Overview
These scripts were tested using Ubuntu 20.04 LTS.

If you are using Windows, you will need to either use specific websites or develop an alternative method to obtain the expected outputs, which some suggestions will be given below.
Another alternative is running a virtual machine with a Linux distro (such as a Long-Term Service (LTS) version of [Ubuntu](https://ubuntu.com/)) or emulation through [Cygwin](https://www.cygwin.com/).


## 01_obtain-kmer-count-via-EMBOSS-compseq.bash
```
USAGE:  01_*.bash  FASTA_FILE  [OUTPUT_DIRECTORY = oligonucleotide_count]
```

This BASH script that goes through the goes through the given `FASTA_FILES` and runs EMBOSS' _compseq_ to obtain the nucleotide motif counts for the given word size (default is using just the forward sense for _k_=1 and both sense and antisense directions for _k_>1).
The output are the _k_-mer files that starts with `[OUTPUT_DIRECTORY]/k#.___`.

For the second script to run successfully, it is critcal that the output starts with `k[_k_-mer]` in the OUTPUT_DIRECTORY. If you are doing an alternative method, you can rename it to whatever you want.

### Alternatives for 01*.bash

You can run compseq on your own computer by downloading the EMBOSS suite (see links in the [Requirements](https://github.com/MoezV/oNTORRA#requirements) section).
Alternatively, you can use a website tat allows you to do so, such as from [Wageningen University & Research (WUR)](https://emboss.bioinformatics.nl/cgi-bin/emboss/compseq).

It is ___critical___ that you account for both forward and reverse directions (i.e. the sense and antisense) for the nucleotide occurrences, which is easily done by _compseq_, to account for directional bias.
The oNTORRA script is set up such that this is not needed for _k_=1, which is accounted for easily, but requires _k_>1 to have both directions accounted (as a simple check of AA==TT, AC==GT, etc. may not be equal, especially at higher _k_ values, due to non-symmetric size - but should be close enough).

## Example
The following is the command line invoked through the BASH script for invoking _compseq_
```bash
compseq -sequence <SEQUENCE_FILE> -word <KMER> -reverse -calcfreq -outfile <OUTPUT_FILE>

OPTIONS
-------
-sequence <SEQUENCE_FILE> = A file containing no FASTA header, only [a single] sequence
-word KMER = the word size (k)
-reverse = Accounts for both the sense and antisense direction counts
-calcfreq = Calculcates the observed and expected frequencies, though this is not required, the following script (02*.bash) assumes it is enabled
-outfile = Output file
```

The following is a portion of the expected file (after having the first ~13 lines skipped which contains some informative details) for _Escherichia coli_ K-12 ([GCF_000005845.2](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000005845.2))
```bash
# Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
#
AA	677590		0.0729902	0.0605389	1.2056751
AC	512472		0.0552036	0.0624844	0.8834791
...
```

## 02_compile-results-from-EMBOSS.bash

This BASH script simply parses the _compseq_ files in to a tabular form in a tab-separated value (TSV) file that has the _k_-mer motif listed in the columns and the organisms listed in the rows.

To perform the alternative for this script, you will have to either make your own (such as using `python`, or a shell command like `tail -n +17 tail -n+17 <KMER_FILE>` to ignore the first 16 lines, etc.)
Try not to do this manually, it will take too long as each _k_-mer has (4^_k_)+1 combinations (e.g. _k_=1-4 have, respectively, 5/17/65/257 combinations, including the 'Other' field).


Example of the tabular file (in a TSV format):

Organism|AA|AC|...
--------|--|--|---
Bacillus_subtilis_168|832142|389101|...
Escherichia_coli_K-12_MG1655|677590|512472|...
...|...|...|...


## 03_run-oNTORRA.R 

The previous two scripts were setting up preparations for this R script, which is the main oNTORRA script.
