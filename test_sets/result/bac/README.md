# oNTORRA bac Result Directory
To view the SVG files, you can open the file directly in most web browsers to render the output image. The SVG file were selected as the test set result image output as it is a vector file which leads to a crisp display and, since it is rendered through the web browser, you can select the data to view the text value in the heatmap which would not be possible with a bitmap-based image file like PNG.

Note that 'Delta' refers to the data calculated though the absoluite difference of the averaged relative abundance values (for _k_-mer _k_#) between two organisms.

```
*_info.txt.gz = Run information (e.g. runtime, session info, etc.)
*_delta_relAbund_k*.tsv.gz = The delta/absoluite difference of the averaged relative abundance values (for _k_-mer _k_#) between two organisms

debug/freq_nt/*_nt-freq_k*.tsv.gz = Contains the calculated odds ratio frequencies (e.g. f*(XY) = f*(XY)/[f*(X)f*(Y)])
debug/rel_abund/*_relAbund_motif_k*.tsv.gz = Contains the calculated odds ratio frequencies for each nucleotide motif combination (e.g. k=2 -> XY, k=3 -> XYZ, k=4 -> XYZW)

Heatmap_Delta/* = contains the SVG for the non-/clustered heatmaps for the Delta values
Heatmap_Motif/* = contains the SVG for the nucleotide motif heatmap for the specified kmer
NJ_Tree/* = Contains the distance tree files. Can view the bootstrapped (1000-replicates) PHYLIP tree file (.phb) in most tree viewing programs like TreeView and MEGA.

(Note that Dim (dimension) = PC (principal component) in the following text)
PCA/*.svg = The SVG of the Dim1 vs Dim2 PCA plots (*PCA-plot*) and the Scree plot (*PCA-Scree-plot*)
PCA/k*/eigenvalue = Contains the eigenvalue and the relative and accumalative variance % contributed for each dimension
PCA/k*/ind_contrib = The nucleotide motif % contribution for each dimension (among all the organism)
PCA/k*/var_contrib = The organism contributed % for each dimension (among all of the k-sized motif)
PCA/k*/summary = The overall standard deviation, and both proportion and cumulative variance, contributed for each principal component for the plots analyzed
```
