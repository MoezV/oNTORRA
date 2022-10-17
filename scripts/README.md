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

__Ensure that the naming convension of the files are identical with k# where only the # changes__. The naming the of _k_=1 file is the only thing that can be different if needed (e.g. if only obtaining the counts/occurrences in only the sense direction), which can be accounted for automatically, but _k_>1 must include the counts in both (sense/antisense) directions

Example of the tabular file (in a TSV format):

Organism|AA|AC|...
--------|--|--|---
Bacillus_subtilis_168|832142|389101|...
Escherichia_coli_K-12_MG1655|677590|512472|...
...|...|...|...


## 03_run-oNTORRA.R 

## Requirements
The previous two scripts were setting up preparations for this R script, which is the main oNTORRA script.

The script itself will check and install the required libraries first (as listed in the [Requirements](https://github.com/MoezV/oNTORRA#requirements) section). As discussed in the previous section, the expected input is a matrix/tabulated data with the nucleotide motifs in the column, the names of the organisms in the row, and the filename  naming convension nearly are identical with only the _k_-mer value is dfferent (though you can edit the script to account for this if needed). The naming the of _k_=1 file is the only thing that can be different if needed (e.g. if only obtaining the counts/occurrences in only the sense direction), which can be accounted for automatically, but _k_>1 must include the counts in both (sense/antisense) directions.

Ensure all the options are set up before running.

## Options
The options are commented in the R script but this will explain further if needed.
```R
# Report settings 
Report_K_Range           = 2:4      # The order/k-mer range wanted to obtain. k=2 for only dinucleotide; k=seq(2,4,1) for di-, tri-, and tetranucleotide results, etc. [Note that trinucleotide is NOT the same as the codon signature (different denominator and frame shift)]
Extended_Genome          = TRUE     # Set to TRUE if the occurrences obtained are based on an "extended" genome where the complementary (antisense) nucleotide bases are appended to the normal genome or the antisense values are accounted for during the process (to account for directional biases)
Add_Antisense_Occurrence = FALSE    # Only functions if Extended_Genome is set to FALSE. If set to TRUE, this script will simply add the number of occurrences to the antisense/complementary sequences to help account for directional biases. If both extended genome and this variable are set to zero, the asterisk for the function variable will be removed as it is used to denote that the complementary bases were accounted for

# Will require the a few parallel libraries. To install, run:  install.packages(c("parallel","foreach","doParallel"))
Enable_Multithread       = TRUE     # Highly recommended to decrease processing time
Max_Threads              = -2       # If value < 0, it will subtract it from the maximum number of threads detected (e.g. if you have 8 threads and maxThreads=-2, 6 threads will be used).
Cluster_Timeout_in_Hour  = 6        # The timeout (in hours) that the pthreads should disconnect when not in use.

# INPUT FILE DEFINITIONS
inputFileName.dir.prefix <- "data/bac/bac3" # The prefix containing the directory and the filename prefix for the K-mer files read
inputFileName.kN.prefix  <- "-k"            # The k-mer occurrences files are named as [PREFIX]-k#-oligo-freq.tsv
inputFileName.kN.suffix  <- "-oligo-freq.tsv"
inputFileName.kN.range   <- seq(2,4)        # The k-mer range to limit. Must be sequential and reach the  value max(Report_K_Range)!
inputFileName.k1.prefix  <- "-"             # If k=1 has a different naming convention, you can define the different here to reflect [DIR_PREFIX] [K1_PREFIX] 1 [K1_SUFFIX]
inputFileName.k1.suffix  <- "_nt-freq.tsv"


#  OUTPUT DEFINITIONS
output.debug.files        <- TRUE    #  Will keep the intermediate files produced to verify if the data was processed correctly (recommended). If FALSE, the debug folder will be deleted!
output.gzip.debug_files   <- TRUE    #  gzip the generated debug files
output.gzip.result_files  <- TRUE    #  gzip the result files (does not affect the info file)
output.gzip.info_file     <- TRUE    #  gzip the info file
output.encoding           <- "UTF-8" # Used for fileEncoding for the output file

dataset.id                <- basename(inputFileName.dir.prefix) # "bac3" # Identifier used for the output names. By default, it is using the prefix after the directory name (or the directory name itself))
resultFile.path           <- paste0("result/", basename(dirname(inputFileName.dir.prefix)), "/")   # The result path directory with the default file prefix
resultFile.file.prefix    <- paste0("oNTORRA_", dataset.id, "_")         #  it is ./result/[input directory name]/oNTORRA_[ID]_
resultFile.path.default   <- paste0(resultFile.path, resultFile.file.prefix)  # Specifically for the main results, These values needs to be separated for the debugging file locations

# You can ignore the following, defining debugging constants
DEBUG                        <- 2      # Set to 0 to disable, 1 for basic debug outputs, 2+ for further debugging outputs 
DEBUG.Average_RevComp_Motifs <- FALSE  # If set to TRUE, it will force an average of [Motif]+[revcomp(Motif)]. As a reminder for this forced averaging, only half the results are shown, with the the reverse-complement motif removed (as it would have an identical result)
DEBUG.force_complement_k1    <- FALSE  # If set to TRUE, will double the occurrences of the k=1 bases (even if it is already f*(k=1)), by adding the complementary base counts (f*(X)=f(X)+f(inv(X))) even if it may already been accounted for in EMBOSS' compseq. Theoretically, this should not affect calculations because usually this just involves the ratios for the variable involved (f(X)/N_k1 != f*(X)/N_k2 BUT f*(X)/N_k2 == C(f*(X))/C(N_k1))), but on the side of caution, leave this variable disabled.
                                    # If set to FALSE (default), it will check if the ratio of round(N_k2/N_k1) is set to either 1 (= f*(k=1)) or 2 (= N_k1 is likely f(k=1)) which in the latter case will have the f(k=1) converted to f*(k=1) as detailed above
#DEBUG.LOAD_DUMMY_DATA     <- FALSE # "Either set value to FALSE or comment this line to disable dummy data"   # Comment out this line to avoid loading the sample dataset (variable itself is a flag). 
debugFile.path.prefix     <- paste0(resultFile.path,"debug/")  # By default, it is placed in the result folder
debugFile.max_species     <- 0      # A limit on how many examples to display. 0 = unlimited
debugFile.max_columns     <- 0      # A limit on how many combinations/calculations to show. 0 = unlimited/show all

# PCA analysis options
options.Output_PCA        <- TRUE   # Principal Component Analysis of the motifs
PCA.Output.Prefix         <- paste0(resultFile.path,"PCA/",resultFile.file.prefix)
PCA.Image_Output          <- "svg"  # Leave blank to see plot without saving. Options: bmp*, jpeg*, pdf+, png*, svg+, tiff*. * = Bitmap images (recommended: PDF), + = vector images (recommended: SVG). Other alternatives: cairo_pdf, cairo_ps
PCA.Width                 <- 680
PCA.Height                <- 600
PCA.Add_Ellipsis          <- FALSE  # Add ellipses to circle each dataset. Default confidence interval: 95% (search for 'plotPCA' to adjust any fine settings)
PCA.Ellipsis.CI           <- 0.95   # Confidence interval for the ellipses.
PCA.Function              <- "prcomp" # prcomp vs princomp


# Disance (NJ) Tree
options.Output_NJ_Tree    <- TRUE   # Print the tree based on the Delta/Diff. Rel. Abund.
                                    # Note that a *.phb (bootstrapped Newick tree) is generated regardless of the image type selected
Tree.Image_Output         <- ""     # Leave blank to disable secondary plot output. Options: bmp*, jpeg*, pdf+, png*, svg+, tiff*. * = Bitmap images (recommended: PDF), + = vector images (recommended: SVG). Other alternatives: cairo_pdf, cairo_ps
Tree.Root_To              <- ""     # Row name of what to root to (leave blank for unrooted)
Tree.Bootstrap            <- 1000   # The number of bootstraps needed
Tree.Output.Prefix        <- paste0(resultFile.path,"NJ_Tree/",resultFile.file.prefix)
                                    #    Appends _unrooted.phb/_rooted-[root].phb depending on the Tree.Root_To value

# Heatmap of the nucleotide motifs
options.Output_Heatmap_Motifs  <- TRUE   # Heatmap based on the Rel. Abund. motifs
HeatMotif.Image_Output         <- "svg" # Leave blank to disable file output. Options: png, pdf, tiff, bmp, jpeg
HeatMotif.Output.Prefix        <- paste0(resultFile.path,"Heatmap/motif/",resultFile.file.prefix)
HeatMotif.Avg_RevComp_Motifs   <- TRUE  # If TRUE, will average the reverse complementary motifs (uses less space)
HeatMotif.Grouping_File        <- "data/bac/_bac-classes.tsv" # Set to empty string "" if not used. Must have headers (col1 header does not matter but col2+ header must reflect the group term you wish to have). At least two columns: col1 = rowname, col2+ = whatever data you wish to reflect (e.g. group name, family name, etc.). The graph reflects the inverse order (i.e. last column would be the first/left-most legend displayed, and conversely, col2 would be the last of the entries)
HeatMotif.Width                <- 1280 #1064
HeatMotif.Height               <- 860  #728

# Heatmap of the absolute difference of the averaged odds ratio between any two organisms
options.Output_Heatmap_Delta   <- TRUE     # Heatmap based on the Rel. Abund. absolute difference between organisms
HeatDelta.Image_Output         <- "svg"    # Leave blank to disable file output. Options: png, pdf, tiff, bmp, jpeg
HeatDelta.Output.Prefix        <- paste0(resultFile.path,"Heatmap/delta/",resultFile.file.prefix)
HeatDelta.Grouping_File        <- "data/bac/_bac-classes.tsv" # Set to empty string "" if not used. Must have headers (col1 header does not matter but col2+ header must reflect the group term you wish to have). At least two columns: col1 = Must reflect rowname, col2+ = whatever data you wish to reflect (e.g. group name, family name, etc.)
HeatDelta.Width                <- HeatMotif.Width  # 1064 px
HeatDelta.Height               <- HeatMotif.Height # 728 px
HeatDelta.Value_Multiplier     <- 1000     # As the values are so small, multiply for readability
HeatDelta.Cluster_Data         <- TRUE     # Enable/disable clustering of the heatmap
HeatDelta.Order_Data_Col       <- 0        # Only applicable if HeatDelta.Cluster_Data is set to FALSE. Sorts the entries based on the column order given. Set to a negative value to disable sorting. If value is zero, it will be ordered using the HeatDelta.Grouping_File columns from the rightmost columns to column 2 in order.

```
