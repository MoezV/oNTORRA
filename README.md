# oNTORRA: OligoNucleoTide Odds Ratio (Relative Abundance), version 0.1.10-R
An R script that is used to calculate the relative abundance values for the specified oligonucleotide kmers and plots the data in several ways:
- Heatmap of the _k_-mer motifs
- Heatmap of the absolute difference (delta) of the averaged relative abundance between several organisms
- Principal Component Analysis (PCA) using the kmer motif
- Distance tree using the Neighbor-Joining algorithm.

# Requirements
1. The oligonucleotide data must first be obtained. It is recommended to use [EMBOSS' compseq](http://emboss.open-bio.org/) ([alt link](https://emboss.bioinformatics.nl/cgi-bin/emboss/compseq)) and generate the occurrences for the word size (_k_-mer, _k_) from _k_=1 to the highest _k_-mer of interest (at least 2) ___in both the sense and anti-sense directions*___. A word size of two (_k_=2) is required to generate the genomic signature (may also be referred simply as either the compositional or dinucleotide bias).
* ___*Important note:___ Ensure that the option to account for the occurrences of the _k_-mer in both the sense and anti-sense directions to account for directional bias. 
    * For the _k_=1 counts, it is not necessary to have this performed as the script will easily account for the anti-sense when the ratio `round(N_k2 / N_k=1) ~ 2`, where N is the total number of occurrences for the given _k_-mer values, or if `DEBUG.force_complement_k1 = TRUE`, but it is absolutely required for _k_ > 1 as the script does NOT account for this nor has an option implemented to do so. It is easily possible to account for this afterwards by programming in to add the occurrences of the  reverse complementary sequence but it is best to take care of this at step 1.
     * To do this:
         * On the [WUR _compseq_ page](https://emboss.bioinformatics.nl/cgi-bin/emboss/compseq), set "Count words in the forward and reverse sense?" to 'Yes'
         * In the _compseq_ software of the [EMBOSS package](ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz), include the flag ``-reverse``

2. Parse the compseq (or an alternative program/script) results such that a tab-separated value (TSV) containing the identifier (e.g. full name of the organism/genome) in the ``row`` and columns with the _k_-mer motifs, with the occurrence values populated. E.g.
    * The filenames of these occurrence files MUST have the same naming convension [with the _k_ # only differing], e.g. `"[set_id]-k#-oligo-freq.tsv"`
    * The _k_=1 filename can differ from _k_>1 and is encouraged to do so if the _k_=1 does not have the anti-sense occurrences accounted for, e.g. `"[set_id]-k1_nt-freq_only-sense-dir.tsv"`

Organism|AA|AC|...
--------|--|--|---
Bacillus_subtilis_168|832142|389101|...
Escherichia_coli_K-12_MG1655|677590|512472|...
...|...|...|...

3. The following R libraries must be installed but you do not need to do this manually as the R script will check if they are present, and if not, will invoke `install.packages()`. The R script was tested using [Rstudio 2022.07.2 Build 576](https://www.rstudio.com/products/rstudio/download/) using [R 4.2.1 (Windows)](https://cran.r-project.org/bin/windows/base/) (see the session info at the bottom of the ```test_set/result/*/oNTORRA_*_info.txt.gz``` file for further details in regards to the test runs.)

R package|Version tested
---------|--------------
parallel (integrated with R-core)|4.1.0
[foreach](https://cran.r-project.org/web/packages/foreach/index.html)|1.5.2
[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)|1.0.17
[factoextra](https://cran.r-project.org/web/packages/factoextra/index.html)|1.0.7
[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)|3.3.6
[ape (Analyses of Phylogenetics and Evolution)](https://cran.r-project.org/web/packages/ape/index.html)|5.6-2
[pheatmap (Pretty Heatmaps)](https://cran.r-project.org/web/packages/pheatmap/index.html)|1.0.12
[RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)|1.1-3
[svglite](https://cran.r-project.org/web/packages/svglite/index.html)|2.1.0


# Citation
If you use this script in your project, please cite this repository and note any major modifications (if applicable) that were made.

There are several ways to cite a GitHub repository, which depends on the journals you are submitting to, but most follow a similar trend as the following:
```
Valliani, Moez. 2022. oNTORRA: OligoNucleoTide Odds Ratio (Relative Abundance), version 0.1.10-R. https://github.com/MoezV/oNTORRA. Date accessed: [date you obtained this repo]
```

It is highly recommended you use a [Referencing Software](https://en.wikipedia.org/wiki/Reference_software) if you are writing a manuscript as it makes the citation process much easier. 
