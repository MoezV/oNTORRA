# compseq results
This directory contains *.tar files which contains the compseq results.

File tree:
```
[set_id]/ = dataset id (e.g. bac)
- 1_nt-counts = containing nucleotide occurrences
-- [genus]/
--- [strain]/
---- *_assembly_report.txt = The genome information that compseq was performed on
---- oligonucleotide_count/
----- k1.*.tsv = The nucleotide counts ONLY IN THE SENSE DIRECTION (f(X))
----- k#-reverse.*.tsv = The nucleotide motif (kmer k#) counts IN BOTH THE SENSE AND ANTI-SENSE DIRECTION (f*(XY...))```
