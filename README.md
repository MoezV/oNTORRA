# oNTORRA: oligoNucleoTide Odds Ratio (Relative Abundance) compositional bias analysis, version 0.1.10-R
[![DOI](https://zenodo.org/badge/552461524.svg)](https://zenodo.org/badge/latestdoi/552461524)

An R script that is used to calculate the relative abundance values for the specified oligonucleotide kmers and plots the data in several ways:
- Heatmap of the _k_-mer motifs
- Heatmap of the absolute difference (delta) of the averaged relative abundance between several organisms
- Principal Component Analysis (PCA) using the kmer motif
- Distance tree using the Neighbor-Joining algorithm.

Please note that this script does __NOT__ calculate the codon signature, the _k_=3 rel. abundance and the codon signature values are different, primarily due to the codon signature having a frame shift of three (for each f(X1,Y2,Z3)) whereas this script only has a frame shift of one. For more information, please see [Karlin et al., 1998](https://github.com/MoezV/oNTORRA/edit/main/README.md#recommended-literature).

# Running oNTORRA
oNTORRA is an R script that parses the counts/occurrences of each oligo-/short nucleotide motifs (with a word size of _k_ (_k_-mer), typically _k_=1-4) in to various display style, primarily a heatmap, PCA analysis, and a distance tree.

To run oNTORRA, please see the [scripts](https://github.com/MoezV/oNTORRA/tree/main/scripts) directory, but ensure you look at the Requirements section below first.

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

# Recommended literature

For more information in regards to the compositional/nucleotide bias analyses, here are several literatures I recommend:
1. Burge, C., Campbell, A.M., and Karlin, S. 1992. _Over- and under-representation of short oligonucleotides in DNA sequences_. PNAS 89: 1358–1362. National Academy of Sciences. [doi:10.1073/pnas.89.4.1358](https://doi.org/10.1073/pnas.89.4.1358).
2. Karlin, S., Campbell, A.M., and Mrázek, J. 1998. _Comparative DNA analysis across diverse genomes_. Annu Rev Genet 32: 185–225. [doi:10.1146/annurev.genet.32.1.185](https://doi.org/10.1146/annurev.genet.32.1.185).
3. Elhai, J. 2001. _Determination of bias in the relative abundance of oligonucleotides in DNA sequences_. J Comput Biol 8: 151–175. [doi:10.1089/106652701300312922](https://doi.org/10.1089/106652701300312922).
4. Phillips, G.J., Arnold, J., and Ivarie, R. 1987. _Mono-through hexanucleotide composition of the Escherichia coli genome: a Markov chain analysis_. Nucleic Acids Res 15: 2611–2626. Oxford Academic. [doi:10.1093/nar/15.6.2611](https://doi.org/10.1093/nar/15.6.2611).
5. Duhagon, M.A., Smircich, P., Forteza, D., Naya, H., Williams, N., and Garat, B. 2011. _Comparative genomic analysis of dinucleotide repeats in Tritryps_. Gene 487: 29–37. [doi:10.1016/j.gene.2011.07.022](https://doi.org/10.1016/j.gene.2011.07.022).
6. Karlin, S., Ladunga, I., and Blaisdell, B.E. 1994. _Heterogeneity of genomes: measures and values_. Proc Natl Acad Sci USA 91: 12837–12841. [doi:10.1073/pnas.91.26.12837](https://doi.org/10.1073/pnas.91.26.12837)
7. Karlin, S., and Mrázek, J. 1997. _Compositional differences within and between eukaryotic genomes_. Proceedings of the National Academy of Sciences 94: 10227–10232. Proceedings of the National Academy of Sciences. [doi:10.1073/pnas.94.19.10227](https://doi.org/10.1073/pnas.94.19.10227).


# Citation
If you use this script in your project, please cite this repository and note any major modifications (if applicable) that were made.

There are several ways to cite a GitHub repository, which depends on the journals you are submitting to. It is highly recommended you use a [Referencing Software](https://en.wikipedia.org/wiki/Reference_software) if you are writing a manuscript as it makes the citation process much easier. 

A few common citation format:

RIS citation
```
TY  - COMP
TI  - oNTORRA: oligoNucleoTide Odds Ratio (Relative Abundance) compositional bias analysis
AU  - Valliani, Moez
AB  - OligoNucleoTide Odds Ratio (Relative Abundance) compositional bias analysis.

Calculates the relative abundance values and plots the data using heatmaps, PCA analysis, and a distance tree.
DA  - 2022/10/17
PY  - 2022
ET  - 0.1.10-R
LA  - R
ST  - oNTORRA
DO  - 10.5281/zenodo.8056684
UR  - https://github.com/MoezV/oNTORRA
ER  - 
```

BibLaTex citation
```
@software{valliani_ontorra,
	title = {{oNTORRA}: {oligoNucleoTide} {Odds} {Ratio} ({Relative} {Abundance}) compositional bias analysis},
	copyright = {GNU GPL v3.0},
	shorttitle = {{oNTORRA}},
	url = {https://github.com/MoezV/oNTORRA},
	version = {0.1.10-R},
	doi = {10.5281/zenodo.8056684},
	abstract = {OligoNucleoTide Odds Ratio (Relative Abundance) compositional bias analysis.

Calculates the relative abundance values and plots the data using heatmaps, PCA analysis, and a distance tree.},
	author = {Valliani, Moez},
	month = oct,
	year = {2022},
}
```

General style citation
```
Valliani, Moez. 2022. oNTORRA: oligoNucleoTide Odds Ratio (Relative Abundance) compositional bias analysis, version 0.1.10-R. https://github.com/MoezV/oNTORRA. DOI: 10.5281/zenodo.8056684. Date accessed: [date you obtained this repo]
```
