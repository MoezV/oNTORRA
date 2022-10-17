# OligoNucleoTide Odds Ratio (Relative Abundance) [oNT-OR(RA) = oNTORRA]
# Version 0.1.10-R (2022-10-10)
# Moez Valliani  2022/09/22

cat("::: Starting oNTORRA run, loading libraries :::\n")
#####   REQUIRED LIBRARIES  ######

# Parallel tutorial: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
for(reqLib in c(
  "parallel",    # version tested: 4.1.0 (should be included in R since 2011)
  "foreach",     # Parallel foreach process. Version used: 1.5.2
  "doParallel",  # To combine with foreach. Version used: 1.0.17
  "factoextra",  # For PCA, Version used: 1.0.7
  "ggplot2",     # For PCA, Version used: 3.3.6
  "ape",         # For producing NJ tree. Version used: 5.6-2. APE = Analyses of Phylogenetics and Evolution
  "pheatmap",    # Version 1.0.12
  "RColorBrewer",# Version 1.1-3
  "svglite"      # Version 2.1.0
)){
  if(library(reqLib, character.only=TRUE, logical.return=TRUE)==FALSE){
    print(paste("Installing",reqLib))
    install.packages(reqLib)
    library(reqLib, character.only=TRUE)
  }
};rm(reqLib)

#####   GLOBAL VARIABLE DEFINITIONS    #####
cat("::: Setting variables :::\n")

Report_K_Range           = 2:4      # The order/k-mer range wanted to obtain. k=2 for only dinucleotide; k=seq(2,4,1) for di-, tri-, and tetranucleotide results, etc. [Note that trinucleotide is NOT the same as the codon signature (different denominator and frame shift)]
Extended_Genome          = TRUE     # Set to TRUE if the occurrences obtained are based on an "extended" genome where the complementary (antisense) nucleotide bases are appended to the normal genome or the antisense values are accounted for during the process (to account for directional biases)
Add_Antisense_Occurrence = FALSE    # Only functions if Extended_Genome is set to FALSE. If set to TRUE, this script will simply add the number of occurrences to the antisense/complementary sequences to help account for directional biases. If both extended genome and this variable are set to zero, the asterisk for the function variable will be removed as it is used to denote that the complementary bases were accounted for

# Will require the a few parallel libraries. To install, run:  install.packages(c("parallel","foreach","doParallel"))
Enable_Multithread       = TRUE     # Highly recommended to decrease processing time
Max_Threads              = -2       # If value < 0, it will subtract it from the maximum number of threads detected (e.g. if you have 8 threads and maxThreads=-2, 6 threads will be used).
Cluster_Timeout_in_Hour  = 6        # The timeout (in hours) that the pthreads should disconnect when not in use.


### INPUT DEFINITIONS
# Expected that the filenames contains the k-mer value and that the filenames are the same (with k1 being the only one that can have another name defined separately)
# E.g.
#      k=1  --> setFolder/setName-[#]-nt-freq.tsv (in this example, k=1 filename has a different suffix)
#      k=2+ --> setFolder/setName-k[#]-oligo-freq.tsv

setwd("X:/path/to/working/directory") # Set your current working directory if needed

inputFileName.dir.prefix <- "data/bac/bac3" # Dataset file prefix
inputFileName.kN.prefix  <- "-k"            # The k-mer occurrences files are named as [SET]-k#-oligo-freq.tsv
inputFileName.kN.suffix  <- "-oligo-freq.tsv"
inputFileName.kN.range   <- seq(2,4)        # The k-mer range to limit. Must be sequential and reach the  value max(Report_K_Range)!
inputFileName.k1.prefix  <- "-"             #  The file containing the single nucleotide occurrences may have a different name associated than the other k-mer files
inputFileName.k1.suffix  <- "_nt-freq.tsv"

# K-MER OCCURRENCES FILE FORMAT
# Header:  Organism	[A* T* C* G* oligomers] Other
# (All files MUST have the same organism names. The others field (which contains any non-ACGT characters) are ignored)


###  OUTPUT DEFINITIONS
output.debug.files        <- TRUE   #  Will keep the intermediate files produced to verify if the data was processed correctly (recommended). If FALSE, the debug folder will be deleted!
output.gzip.debug_files   <- TRUE   #  gzip the generated debug files
output.gzip.result_files  <- TRUE  #  gzip the result files (does not affect the info file)
output.gzip.info_file     <- TRUE   #  gzip the info file
output.encoding           <- "UTF-8" # Used for fileEncoding for the output file

dataset.id                <- basename(inputFileName.dir.prefix) # "bac3" # Identifier used for the output names. By default, it is using the prefix after the directory name (or the directory name itself))
resultFile.path           <- paste0("result/", basename(dirname(inputFileName.dir.prefix)), "/")   # The result path directory with the default file prefix
resultFile.file.prefix    <- paste0("oNTORRA_", dataset.id, "_")         #  it is ./result/[input directory name]/oNTORRA_[ID]_
resultFile.path.default   <- paste0(resultFile.path, resultFile.file.prefix)  # Specifically for the main results, These values needs to be separated for the debugging file locations

# You can ignore the following, defining constants
DEBUG                        <- 2      # Set to 0 to disable, 1 for basic debug outputs, 2+ for further debugging outputs 
DEBUG.Average_RevComp_Motifs <- FALSE  # If set to TRUE, it will force an average of [Motif]+[revcomp(Motif)]. As a reminder for this forced averaging, only half the results are shown, with the the reverse-complement motif removed (as it would have an identical result)
DEBUG.force_complement_k1    <- FALSE  # If set to TRUE, will double the occurrences of the k=1 bases (even if it is already f*(k=1)), by adding the complementary base counts (f*(X)=f(X)+f(inv(X))) even if it may already been accounted for in EMBOSS' compseq. Theoretically, this should not affect calculations because usually this just involves the ratios for the variable involved (f(X)/N_k1 != f*(X)/N_k2 BUT f*(X)/N_k2 == C(f*(X))/C(N_k1))), but on the side of caution, leave this variable disabled.
                                    # If set to FALSE (default), it will check if the ratio of round(N_k2/N_k1) is set to either 1 (= f*(k=1)) or 2 (= N_k1 is likely f(k=1)) which in the latter case will have the f(k=1) converted to f*(k=1) as detailed above
#DEBUG.LOAD_DUMMY_DATA     <- FALSE # "Either set value to FALSE or comment this line to disable dummy data"   # Comment out this line to avoid loading the sample dataset (variable itself is a flag). 
debugFile.path.prefix     <- paste0(resultFile.path,"debug/")  # By default, it is placed in the result folder
debugFile.max_species     <- 0      # A limit on how many examples to display. 0 = unlimited
debugFile.max_columns     <- 0      # A limit on how many combinations/calculations to show. 0 = unlimited/show all

options.Output_PCA        <- TRUE   # Principal Component Analysis of the motifs
PCA.Output.Prefix         <- paste0(resultFile.path,"PCA/",resultFile.file.prefix)
PCA.Image_Output          <- "svg"  # Leave blank to see plot without saving. Options: bmp*, jpeg*, pdf+, png*, svg+, tiff*. * = Bitmap images (recommended: PDF), + = vector images (recommended: SVG). Other alternatives: cairo_pdf, cairo_ps
PCA.Width                 <- 680
PCA.Height                <- 600
PCA.Add_Ellipsis          <- FALSE  # Add ellipses to circle each dataset. Default confidence interval: 95% (search for 'plotPCA' to adjust any fine settings)
PCA.Ellipsis.CI           <- 0.95   # Confidence interval for the ellipses.
PCA.Function              <- "prcomp" # prcomp vs princomp
# Depreciated. Set up so it does both regular PCA1v2 and Screen plots. # PCA.Scree_Plot            <- FALSE  # Scree plot - visualized eigenvalues

options.Output_NJ_Tree    <- TRUE   # Print the tree based on the Delta/Diff. Rel. Abund.
                                    # Note that a *.phb (bootstrapped PHYLIP tree) is generated regardless of the image type selected
Tree.Image_Output         <- ""     # Leave blank to disable secondary plot output. Options: bmp*, jpeg*, pdf+, png*, svg+, tiff*. * = Bitmap images (recommended: PDF), + = vector images (recommended: SVG). Other alternatives: cairo_pdf, cairo_ps
Tree.Root_To              <- ""     # Row name of what to root to (leave blank for unrooted)
Tree.Bootstrap            <- 1000   # The number of bootstraps needed
Tree.Output.Prefix        <- paste0(resultFile.path,"NJ_Tree/",resultFile.file.prefix)
                                    #    Appends _unrooted.phb/_rooted-[root].phb depending on the Tree.Root_To value

options.Output_Heatmap_Motifs  <- TRUE   # Heatmap based on the Rel. Abund. motifs
HeatMotif.Image_Output         <- "svg" # Leave blank to disable file output. Options: png, pdf, tiff, bmp, jpeg
HeatMotif.Output.Prefix        <- paste0(resultFile.path,"Heatmap/motif/",resultFile.file.prefix)
HeatMotif.Avg_RevComp_Motifs   <- TRUE  # If TRUE, will average the reverse complementary motifs (uses less space)
HeatMotif.Grouping_File        <- "data/bac/_bac-classes.tsv" # Set to empty string "" if not used. Must have headers (col1 header does not matter but col2+ header must reflect the group term you wish to have). At least two columns: col1 = rowname, col2+ = whatever data you wish to reflect (e.g. group name, family name, etc.). The graph reflects the inverse order (i.e. last column would be the first/left-most legend displayed, and conversely, col2 would be the last of the entries)
HeatMotif.Width                <- 1280 #1064
HeatMotif.Height               <- 860  #728


options.Output_Heatmap_Delta   <- TRUE     # Heatmap based on the Rel. Abund. absolute difference between organisms
HeatDelta.Image_Output         <- "svg"    # Leave blank to disable file output. Options: png, pdf, tiff, bmp, jpeg
HeatDelta.Output.Prefix        <- paste0(resultFile.path,"Heatmap/delta/",resultFile.file.prefix)
HeatDelta.Grouping_File        <- "data/bac/_bac-classes.tsv" # Set to empty string "" if not used. Must have headers (col1 header does not matter but col2+ header must reflect the group term you wish to have). At least two columns: col1 = Must reflect rowname, col2+ = whatever data you wish to reflect (e.g. group name, family name, etc.)
HeatDelta.Width                <- HeatMotif.Width  # 1064 px
HeatDelta.Height               <- HeatMotif.Height # 728 px
HeatDelta.Value_Multiplier     <- 1000     # As the values are so small, multiply for readability
HeatDelta.Cluster_Data         <- TRUE     # Enable/disable clustering of the heatmap
HeatDelta.Order_Data_Col       <- 0        # Only applicable if HeatDelta.Cluster_Data is set to FALSE. Sorts the entries based on the column order given. Set to a negative value to disable sorting. If value is zero, it will be ordered using the HeatDelta.Grouping_File columns from the rightmost columns to column 2 in order.

# Generating the folders for the outputs of interests 
## Result path
resultOutput.info         <- paste0(resultFile.path.default, "info.txt")


# Other variables
options.RNG_Seed = 4                             # Setting seed for reproducability 
Pixels_Per_Inch = 76                             # Roughly 76 pixels per inch (required for outputting vector files like PDF or SVG formats)
knt_letters = c(LETTERS[24:26], LETTERS[23:1])   #  Letters used for general pattern expression (XYZ WVU...A)
invNT = c("A"="T","T"="A","C"="G","G"="C")
Regex_Remove_Filename_Chars="[^0-9A-Za-z._:/ -]" #  Regular expression script to ensure special characters are removed using gsub
                                                ## gsub(Regex_Remove_Filename_Chars, "", [string], ignore.case=TRUE)
# Ignored/future implantation:
#  codonFreq = NOT IMPLEMENTED. Get the codon-frequency rather than mononucleotide frequency and treating each position/word size a block of 3 instead of 1 (e.g. rho*(x1,y2,z3))
#  parseFreq = NOT IMPLEMENTED. Use the frequencies calculated from EMBOSS rather than calculating it here with the occurances observed.
#  Use_HOMM  = NOT IMPLEMENTED. If set to TRUE, any value for Report_K_Range above the given inputFileName.kN.range below will be calculated using the Highest Order Markov Chain. E.g. If the report goes from k=2-6 but the observed frequencies were only given for k=2-4, then k=5-6 will be calculated with the fourth-order Markov model.

### END OF GLOBAL OPTIONS ###
#####  MAIN SCRIPT #####

######  Sanitizing input values (safety check 1)  ######

# kmer range must have a value >= 2
if(max(inputFileName.kN.range) < max(Report_K_Range)){ stop("The maximum value for the file input range [max(inputFileName.kN.range) = ", max(inputFileName.kN.range), "] is less than the maximum value of the Report_K_Range [max(Report_K_Range) = ", max(Report_K_Range),"]") }
if(min(Report_K_Range) < 2){ stop("ERROR: The lowest value (", min(Report_K_Range),") is below the lowest threshold (min(k) >= 2)!") }
if(max(Report_K_Range) > 4){ warning("oNTORRA has only been tested for k=2-4. Theoretically, higher k-sizes should be fine based on debugging tests but you should verify the formula output in case!") }
if(Extended_Genome==TRUE && Add_Antisense_Occurrence==TRUE){ stop("Please double check your configurations. Cannot have the genome counts defined as being extended (antisense already accounted for) and wanting to account for the antisense occurrences...") }
if(Extended_Genome==FALSE && Add_Antisense_Occurrence==FALSE){ stop("Please double check your configurations. Frequencies will NOT be accounting for directional bias (HIGHLY RECOMMENDED YOU TURN ON ANTISENSE CALCULATION)! If you do want to go ahead with this, please search for this message and change stop() to warning().") }
if(exists("DEBUG.LOAD_DUMMY_DATA") && get("DEBUG.LOAD_DUMMY_DATA")!=FALSE){
  warning("Warning: dummy data flag present! Dummy data will be used! (To disable, comment out the variable DEBUG.LOAD_DUMMY_DATA)") 
  # Update dataset names for dummy data
  dataset.id  = "dummyData"
  resultFile.path = "result/dummyData/"
  resultFile.file.prefix = "oNTORRA_dummyData_"
  resultFile.path.default = paste0(resultFile.path, resultFile.file.prefix)
  debugFile.path.prefix = paste0(resultFile.path,"debug/")
  resultOutput.info = paste0(resultFile.path.default, "info.txt");
}else{
  resultFile.path = gsub(Regex_Remove_Filename_Chars, "", resultFile.path, ignore.case = TRUE)
  resultFile.file.prefix = gsub(Regex_Remove_Filename_Chars, "", resultFile.file.prefix, ignore.case = TRUE)
  resultFile.file.prefix = gsub(Regex_Remove_Filename_Chars, "", resultFile.file.prefix, ignore.case = TRUE)
  resultFile.path.default = gsub(Regex_Remove_Filename_Chars, "", resultFile.path.default, ignore.case = TRUE)
  resultOutput.info = gsub(Regex_Remove_Filename_Chars, "", resultOutput.info, ignore.case = TRUE)
}
set.seed(options.RNG_Seed)

# Result check path
checkPath = resultFile.path; if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)

## Default debug paths
if(output.debug.files==TRUE){
  # Nucleotide frequency debugging options
  debugOutput.freq.path <- gsub(Regex_Remove_Filename_Chars, "",
                                paste0(debugFile.path.prefix,"freq_nt/"), ignore.case = TRUE);  checkPath = debugOutput.freq.path; if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
  debugOutput.freq.omit_zero_columns <- TRUE  #  If TRUE, all zero-only columns would be omitted. If you want to see all 4^kmer possibilities, set to FALSE.
  
  # Relative abundance debugging options
  debugOutput.relAbund.motif.path <- gsub(Regex_Remove_Filename_Chars, "", 
                                          paste0(debugFile.path.prefix,"rel_abund/"), ignore.case = TRUE); checkPath = debugOutput.relAbund.motif.path; if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
  debugOutput.relAbund.motif.omit_zero_columns <- TRUE  #  If TRUE, all zero-only columns would be omitted. If you want to see all 4^kmer possibilities, set to FALSE.
  
}


###### Important functions #######
cat("::: Setting functions :::\n")

iterN <- function(motif,spos=2, replChar="N"){
  # Iterative function to obtain combinations of N-replaced *inner* positions, e.g. XYZW -> XNZW, XYNW, XNNW
  # Iterative replacement function taken from the snippet posted by moodymudskipper on 2022/09/23 on StackOverflow (https://stackoverflow.com/questions/52200936/create-all-combinations-of-letter-substitution-in-string)
  #
  # Returns:
  #    Iterative replacement of inner-positioned characters, leaving the two outer characters alone
  #
  # ARGUMENTS:
  #    motif    = E.g. "ATG" --> c(ATG, ANG).  "ATCG" --> c(ATCG, ATNG, ANCG, ANNG)
  #    spos     = starting position (increases in each iteration)
  #    replChar = The character replacing the inner positions (default: "N")

  if(length(motif)>1){  # Safety check in case a vector of motifs were submitted, run each iteration separately
    warning("VECTOR MOTIFS IN iterN() [treating each separately]: ",paste(motif, collapse=", "))
    allMotif=c()
    for(Xtra in motif){
      allMotif=c(allMotif, iterN(Xtra, spos, replChar))
    }
    allMotif
  }else{
    if((nchar(motif) > 2) && (spos < nchar(motif)))
      c(iterN(motif,spos+1,replChar), iterN(`substr<-`(motif, spos, spos, replChar),spos+1,replChar))
    else 
      # if(0&&omitInit==TRUE){motif[-1]}
      motif
  }
};if(DEBUG>3){ cat("## Debugging N-replacement iteration. Testing XYZW, expecting: XYNW, XNZW, XNNW  ##\n"); iterN("XYZW")[-1] }

permATCG <- function(klength=2, Derive=FALSE , replChar=".", groupFUN=c){
  # To generate permutations of ACGT of length/word size 'klength'
  #
  # Expectations:
  #    If Derive == FALSE: return 4^klength combinations (only permutations of klength ACGT combinations)
  #    If Derive == TRUE: returns returns 4^(klength!)-ish combinations
  #
  # Returns:
  #    Permutation ACTG for word size 'klength' (and its immediate derivatives if TRUE)
  #
  # Arguments:
  #    klength   = Kmer/word size
  #    Derive    = If lower length combinations with replacements should be taken in to account (klength to 1)
  #    replChar  = Iterative replacement characters for the inner positions of the derived permutations (e.g. XYZW --> X.ZW, XY.W, X..W, X.Z, etc.. If "", iterative replacement is skipped
  #    groupFUN  = Grouping function. Values:
  #                   c    <- vector, everything together rather than split apart
  #                   list <- list, separates everything by kmer-products (useful for debugging)

  ACGT <- c("A","C","G","T")
  outResult = groupFUN()
  for(ksize in klength:1){
    curATCG = apply(expand.grid(rep.int(list(ACGT), ksize)),1,paste0,collapse="")  ## For unique combinations (not permutations): unlist(combn(c("A","C","G","T"), klength, simplify=F, FUN=paste0, collapse=""))
    if(Derive==FALSE){
      outResult = groupFUN(curATCG);
      break;
    }else{
      if(ksize>2 && replChar != ""){
        outResult = groupFUN(outResult, unique(iterN(curATCG, replChar=replChar)))
      }else{
        outResult = groupFUN(outResult, curATCG)
      }
    }
  }
  unique(outResult, rm.na=T)
};if(DEBUG>3){ print("### Debugging permATCG"); tmp <- permATCG(4, Derive=T, groupFUN=list); print(tmp); print(length(tmp)); rm(tmp) }

subsetNT<-function(motif, kmer, iterMotif=FALSE, replChar="."){
  if(kmer==1){return(strsplit(motif,"")[[1]])}
  kcombin = c()
  for(epos in 1:nchar(motif)){ # Go from the current iteration to the end of the kmer/oligonucleotide word size.
    for(spos in 1:nchar(motif)){
      if(spos >= epos){ next }
      curSeq = substr(motif, spos, epos)
      if(nchar(curSeq) != kmer){next}

      # If the extraction is much larger than the expected iteration, it means there are components in the middle that needs to be replaced with N
      if(iterMotif){ kcombin=c(kcombin,iterN(curSeq, replChar = replChar)) }
      kcombin=c(kcombin,curSeq)
      
    }
  }
  kcombin
}; if(DEBUG>3){ print(subsetNT("XYZW", 2)) }

derivATCG <- function(motif, iterInner=TRUE, replChar=".", groupFUN=c, splitOddEven=TRUE){
  # Generates derivations of 'motif' with the inner-characters iteratively replaced with ./N
  #
  # Returns:
  #    list(numerator derived motifs, denominator derived motifs)
  #
  # Arguments:
  #    motif      = E.g. "ATC" to obtain c(ATC, A.C, AT, TC, A, T, C)
  #    iterInner = If set to FALSE, inner replacements will not be done (thus end result would be, for example, "ATC" --> c(ATC,AT,TC,A,T,C))
  #    replChar  = Iterative replacement characters for the inner positions of the derived permutations (e.g. XYZW --> X.ZW, XY.W, X..W, X.Z, etc.. If "", iterative replacement is skipped
  #    groupFUN  = Grouping function. Values:
  #                   c    <- vector, everything together rather than split apart
  #                   list <- list, separates everything by kmer-products (useful for debugging)
  
  ACGT <- c("A","C","G","T")
  motif = toupper(motif)  # Ensure sequences are uppercase for the replacement
  outResult = groupFUN()
  firstLoop=TRUE
  for(ksize in nchar(motif):1){
    if(firstLoop==TRUE){
      firstLoop=FALSE
      curATCG = motif
    }else{
      curATCG = subsetNT(motif, ksize) #apply(expand.grid(rep.int(strsplit(motif,""), ksize)),1,paste0,collapse=""); # unlist(permn(c("A","T","C"), paste0, collapse=""))
    }
    if(ksize>2 && iterInner==TRUE && replChar != ""){
      outResult = groupFUN(outResult, (iterN(curATCG, replChar=replChar))) # unique
    }else{
      outResult = groupFUN(outResult, curATCG)
    }
  }
  ### Want to avoid unique() so things like TTTT can output two T.T possibilities (T1.T3, T2.T4) :: 
  ###outResult = unique(outResult)
  if(splitOddEven==TRUE){ 
    numChar = nchar(gsub("[.]","",outResult))  #  Treat N/. positions as not being present
    revOrder = rev(order(numChar))
    revChar = numChar[revOrder]
    rMotif = outResult[revOrder]
    splitFirst = ((revChar%%2)==(nchar(motif)%%2))
    list(rMotif[splitFirst],rMotif[!splitFirst])
  }else{
    outResult
  }
};if(DEBUG>3){ print("### Debugging derivATCG");  tmp <- derivATCG("XYZW", iterInner=T, groupFUN=c); print(tmp); print(length(tmp)); rm(tmp) }

revcomp<-function(motif, Nchar="."){
  # A function to reverse-complement NT strings
  # Nchar = non-ATCG character replacements
  motif = toupper(motif) # Ensuring replacement will occur correctly
  newMotif=""
  invNT = c("A"="T","T"="A","C"="G","G"="C","N"=Nchar,"."=Nchar)
  for(NT in strsplit(motif,"")[[1]]){  # Perform the complementary replacement by each character
    newMotif=paste0(invNT[NT], newMotif) #str_replace_all(NT, c("A"="T","T"="A","C"="G","G"="C","[^ATCG]"=Nchar)), newMotif)
  }
  return(newMotif)
}


######   Main script   ######


######  Parallel parameters ######
{
  # Enabling multi-threading. Usage learned from the 'Multi-threading in R' blog post by Jospeh Crispell on August 27, 2018 at https://josephcrispell.github.io/2018/08/27/multi-threading-R.html
  parallel.enabled = FALSE
  parallel.threads = 1
  if(Enable_Multithread==TRUE){
    cat("::: Setting parallel threads (pthreads) :::\n")
    parallel.enabled = TRUE
    
    # Positive max_thread = definite number of cores. Negative value: subtract from total threads present
    if(Max_Threads>0){
      parallel.threads = min(detectCores()-1, Max_Threads)  # Ensure we do not use all the threads otherwise everything will slow down
    }else{
      parallel.threads = detectCores() + Max_Threads
    }
    parallel.threads = max(1, min(parallel.threads, detectCores()-1))  # Ensure there is at least one thread in use and at least one core available (if possible).
  }
  ######## RESTART CLUSTER ########
  if(exists("pClust")==TRUE){ stopCluster(pClust); rm(pClust) }  # Safety check if a parallel cluster is already set
  pClust <- makeCluster(parallel.threads, timeout = (Cluster_Timeout_in_Hour*60*60))    # Make the clusters available. Timeout: 2 hours without response (4 hr * 60 min/hr * 60 sec/min = 14,400 sec)
}


### LOAD/CREATE DATA IN TO NEW VAR:  occ.k[#]

cat("::: Loading files :::\n")

####### Dummy data #######
# Check if the dummy data should be loaded (dummy seq:  AATCCGGAT)
if(exists("DEBUG.LOAD_DUMMY_DATA") && get("DEBUG.LOAD_DUMMY_DATA")!=FALSE){ # In case DEBUG.LOAD_DUMMY_DATA is defined with FALSE instead of being uncommented
  # NOTE: DUMMY SET FREQUENCY VALUES ARE FOUND AT THE END OF THE SCRIPT!
  ACGT = c("A", "C", "G", "T")
  ACGTX = c("A", "C", "G", "T", "Other")
  occ.k1 = matrix(data = 0, nrow = 1, ncol = 5, dimnames = list("null", ACGTX)); occ.k1
  occ.k2 = cbind(matrix(data = 0, nrow = 1, ncol = 4**2, dimnames = list("null", apply(expand.grid(rep.int(list(ACGT),2)),1,paste0,collapse=""))),0); colnames(occ.k2)[length(occ.k2)]="Other"; occ.k2
  occ.k3 = cbind(matrix(data = 0, nrow = 1, ncol = 4**3, dimnames = list("null", apply(expand.grid(rep.int(list(ACGT),3)),1,paste0,collapse=""))),0); colnames(occ.k3)[length(occ.k3)]="Other"; occ.k3
  occ.k4 = cbind(matrix(data = 0, nrow = 1, ncol = 4**4, dimnames = list("null", apply(expand.grid(rep.int(list(ACGT),4)),1,paste0,collapse=""))),0); colnames(occ.k4)[length(occ.k4)]="Other"; occ.k4
  
  
  warning("WARNING: Loading dummy dataset 1 [AATCCGGAT]! [f(k=1), f*(k>1)]")
  occ.k1 = rbind(occ.k1,matrix(data = c(3,2,2,2,0), nrow = 1, ncol = 5, dimnames = list("dummy_set_1", ACGTX))); occ.k1
  occ.k2 = merge(occ.k2,matrix(data = c(1,4,2,2,2,2,2,1), nrow = 1, ncol = 8, dimnames = list("dummy_set_1", c("AA","AT","CC","CG","GA","GG","TC","TT"))),all=T,no.dups=F,sort=F); occ.k2
  occ.k3 = merge(occ.k3,matrix(data = c(1,2,1,2,2,2,2,2), nrow = 1, ncol = 8, dimnames = list("dummy_set_1", c("AAT", "ATC", "ATT", "CCG", "CGG", "GAT", "GGA", "TCC"))),all=T,no.dups=F,sort=F); occ.k3
  occ.k4 = merge(occ.k4,matrix(data = c(1,2,2,2,1,2,2), nrow = 1, ncol = 7, dimnames = list("dummy_set_1", c("AATC", "ATCC", "CCGG", "CGGA", "GATT", "GGAT", "TCCG"))),all=T,no.dups=F,sort=F); occ.k4

  
  warning("WARNING: Loading dummy dataset 2 [CGATCCGGATA]! (f*[k>0])")
  occ.k1 = rbind(occ.k1,matrix(data = c(5,6,6,5,0), nrow = 1, ncol = 5, dimnames = list("dummy_set_2", c("A","C","G","T","Other")))); occ.k1
  occ.k2 = merge(occ.k2,matrix(data = c(4,2,4,3,2,2,3), nrow = 1, ncol = 7, dimnames = list("dummy_set_2", c("AT","CC","CG","GA","GG","TA","TC"))),all=T,no.dups=F,sort=F); occ.k2
  occ.k3 = merge(occ.k3,matrix(data = c(1,3,2,1,2,3,2,1,2,1), nrow = 1, ncol = 10, dimnames = list("dummy_set_2", c("ATA","ATC","CCG","CGA","CGG","GAT","GGA","TAT","TCC","TCG"))),all=T,no.dups=F,sort=F); occ.k3
  occ.k4 = merge(occ.k4,matrix(data = c(2,1,2,1,2,1,2,2,1,2), nrow = 1, ncol = 10, dimnames = list("dummy_set_2", c("ATCC","ATCG","CCGG","CGAT","CGGA","GATA","GATC","GGAT","TATC","TCCG"))),all=T,no.dups=F,sort=F); occ.k4
  

  
  warning("WARNING: Loading dummy dataset 3 [TGACGNATTAACNNNGGATC]! (f*[k>0])")
  
  # f(k=1) set, f*(k=1) is: 9,7,7,9,8
  occ.k1 = rbind(occ.k1,matrix(data = c(5,3,4,4,4), nrow = 1, ncol = 5, dimnames = list("dummy_set_3", ACGTX))); occ.k1
  occ.k2 = merge(occ.k2,matrix(data = c(2,2,0,4,1,1,2,0,3,0,1,2,2,3,1,2,12), nrow = 1, ncol = 17, dimnames = list("dummy_set_3", c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT","Other"))),all=T,no.dups=F,sort=F); occ.k2
  occ.k3 = merge(occ.k3,matrix(data = c(1,1,1,2,1,1,1,2,1,1,1,2,1,1,1,2,16), nrow = 1, ncol = 17, dimnames = list("dummy_set_3", c("AAC","AAT","ACG","ATC","ATT","CGT","GAC","GAT","GGA","GTC","GTT","TAA","TCA","TCC","TGA","TTA","Other"))),all=T,no.dups=F,sort=F); occ.k3
  occ.k4 = merge(occ.k4,matrix(data = c(1,1,1,1,2,1,1,1,1,1,1,2,20), nrow = 1, ncol = 13, dimnames = list("dummy_set_3", c("ATCC","ATTA","CGTC","GACG","GATC","GGAT","GTCA","GTTA","TAAC","TAAT","TGAC","TTAA","Other"))),all=T,no.dups=F,sort=F); occ.k4

  
  # Cleaning up the occ dummy tables
  occ.k1 = occ.k1[-1,]
  occ.k2 = occ.k2[-1,];occ.k2=apply(occ.k2,2,function(x){ifelse(is.na(x),0,x)});rownames(occ.k2) = rownames(occ.k1); occ.k2
  occ.k3 = occ.k3[-1,];occ.k3=apply(occ.k3,2,function(x){ifelse(is.na(x),0,x)});rownames(occ.k3) = rownames(occ.k1); occ.k3[,1:18]
  occ.k4 = occ.k4[-1,];occ.k4=apply(occ.k4,2,function(x){ifelse(is.na(x),0,x)});rownames(occ.k4) = rownames(occ.k1); occ.k4[,1:18]
  
  # Verifying length. Assuming 'Other' field is present (thus N_k = 4^k + 1)
  for(fi in 1:4){if(dim(get(paste0("occ.k",fi)))[2] != (1+4^fi)){ stop(paste0("ERROR: Dummy set ", fi, "does not match the expected column size! Was 'Other' removed?!")) }}

  cat("::: LOADED DUMMY SET VALUES :::\n")
}else{
  ###### Read Data ######
  # Loading regular data, reading the k-mer (k>1) occurrence files
  for(kmer in inputFileName.kN.range){
    if(!exists(paste0("occ.k", kmer))){
      inputFileName = paste0(inputFileName.dir.prefix, inputFileName.kN.prefix, kmer, inputFileName.kN.suffix)
      if(!file.exists(inputFileName)){ stop(paste0("K1 occurrence file does not exist: '", inputFileName,"'")); }
      
      tmp = paste0("occ.k", kmer)  #  Assigning the new variable names through assign() as they involve variables
      assign(tmp, read.delim(inputFileName, header=TRUE, sep="\t",row.names=1));
      if(DEBUG){
        cat(paste0("--- ", tmp, " ---\n"))
        print(dim(get(tmp)));
        print(get(tmp)[1:3,1:5]);
        cat("\n")
      }
    }
  }
  #  || get("occ.k2")==FALSE
  if(!exists("occ.k2")){ stop("k=2 file was not read or does not exist (e.g. forgot to define in k-range). Need at least k-mer of two for the rel. abundance calculation!") }
  
  # Reading the k1 (single nucleotide [NT] occurrences) if it has a different name defined
  if(!exists("occ.k1") && inputFileName.k1.suffix!=""){
    inputFileName = paste0(inputFileName.dir.prefix, inputFileName.k1.prefix, 1, inputFileName.k1.suffix)
    if(!file.exists(inputFileName)){ stop(paste0("K1 occurrence file does not exist: '", inputFileName,"'")); }
    
    occ.k1 <- read.delim(inputFileName, header=TRUE, sep="\t",row.names=1);
    if(DEBUG){ cat("--- occ.k1 ---\n"); print(dim(occ.k1)); print(occ.k1[1:3,1:5]); cat("\n")}
  }
  #  || get("occ.k1")==FALSE
  if(!exists("occ.k1")){ stop("k=1 file was not read or does not exist (e.g. forgot to define in k-range)") }
  
  cat("::: Loaded dataset :::\n")
}

######  Verifying number of parsed entries (safety check 2)  ######
for(ii in 1:max(Report_K_Range)){
  check_Kmer_file = paste0("occ.k", ii)
  if(exists(check_Kmer_file)==FALSE){
    stop(paste0("ERROR: The k-mer input for k = ", ii, " is not found. Please make sure you have the appropriate range for inputFileName.kN.range (",inputFileName.kN.range,")"))
  }
}; rm(check_Kmer_file, ii)

###### Debug info output ######
# Output the input information to the info file before the variables are deleted
### If gzipped option enabled, set up the persistent gzip connection
if(output.gzip.info_file==TRUE){ resultOutput.info = gzfile(paste0(resultOutput.info,".gz"), open="w", compression = 9, encoding = output.encoding); 
}else{resultOutput.info = file(resultOutput.info, "w", encoding = output.encoding)}

cat(paste0("oNTORRA [ oligoNucleoTide Odds Ratio (Relative Abundance) ]\n",
          "Version: 0.1.10-R (2022-10-10)\n",
          "\U00A9 2022 Moez Valliani\n",
          strrep("-", 57),
          "\n\n",
          
          "___  INPUT INFO  ___\n\t",
          "---  Rel. Abund. ---\n\t",
          "Dataset ID:", strrep("\t",3), dataset.id, "\n\t",
          "Result Directory prefix:", strrep("\t",1),inputFileName.dir.prefix,"\n\t",
          "k1 file prefix/suffix:", strrep("\t",2), inputFileName.k1.prefix, " / ", inputFileName.k1.suffix, "\n\t",
          "k_n file prefix/suffix:", strrep("\t",2), inputFileName.kN.prefix, " / ", inputFileName.kN.suffix, "\n\t",
          "k_n range:", strrep("\t",3), paste(inputFileName.kN.range, collapse=","),"\n\n",
          
          "--- Principal Component Analysis (PCA) ---\n\t",
          "function:", strrep("\t",3),PCA.Function,"\n\t",
          "Enabled:", strrep("\t",3),options.Output_PCA,"\n\t",
          "Output prefix:", strrep("\t",3),PCA.Output.Prefix,"\n\t",
          "Output image extension:", strrep("\t",2),PCA.Image_Output,"\n\t\t",
          "Width:", strrep("\t",3), PCA.Width, "\n\t\t",
          "Height:", strrep("\t",3),PCA.Height, "\n\t",
          "Ellipses added:", strrep("\t",3),PCA.Add_Ellipsis, "\n\t\t",
          "Confidence Interval:", strrep("\t",1),PCA.Ellipsis.CI, "\n\n",
          
          "--- Neighbor-Joining (NJ) Tree ---\n\t",
          "Enabled:", strrep("\t",3),options.Output_NJ_Tree,"\n\t",
          "Output prefix:", strrep("\t",3),Tree.Output.Prefix,"\n\t",
          "Output image extension:", strrep("\t",1),Tree.Image_Output,"\n\t",
          "Rooted to:", strrep("\t",2), Tree.Root_To,"\n\t",
          "Bootstrap:", strrep("\t",3), Tree.Bootstrap,"\n\n",

          "--- Heatmap (Rel. Abund. Motif) ---\n\t",
          "Enabled:", strrep("\t",3),options.Output_Heatmap_Motifs,"\n\t",
          "Output prefix:", strrep("\t",3),HeatMotif.Output.Prefix,"\n\t",
          "Output image extension:", strrep("\t",2),HeatMotif.Image_Output,"\n\t\t",
          "Width:", strrep("\t",3), HeatMotif.Width, "\n\t\t",
          "Height:", strrep("\t",3),HeatMotif.Height, "\n\t",
          "Group file:", strrep("\t",3), HeatMotif.Grouping_File, "\n\n",
          
          "--- Heatmap (Absol. Diff. Rel. Abund.) ---\n\t",
          "Enabled:", strrep("\t",3),options.Output_Heatmap_Delta,"\n\t",
          "Output prefix:", strrep("\t",3),HeatDelta.Image_Output,"\n\t",
          "Output image extension:", strrep("\t",2), HeatDelta.Output.Prefix,"\n\t\t",
          "Width:", strrep("\t",3), HeatDelta.Width, "\n\t\t",
          "Height:", strrep("\t",3),HeatDelta.Height, "\n\t",
          "Average RevComp motifs:", strrep("\t",1), HeatMotif.Avg_RevComp_Motifs, "\n\t",
          "Group file:", strrep("\t",3), HeatDelta.Grouping_File, "\n\n",

          "___  DEBUG INFO  ___\n\t",
          "Debug Directory prefix:", strrep("\t",2),debugFile.path.prefix,"\n\t",
          "Dummy data set loaded:", strrep("\t",2), (exists("DEBUG.LOAD_DUMMY_DATA") && get("DEBUG.LOAD_DUMMY_DATA")!=FALSE), "\n\t",
          "Force complement k=1 values:\t", DEBUG.force_complement_k1, "\n\t",
          "Max species displayed:", strrep("\t",2), debugFile.max_species, ifelse(debugFile.max_species==0," (show all)",""), "\n\t",
          "Max combo displayed:", strrep("\t",3-1), debugFile.max_columns, ifelse(debugFile.max_columns==0," (show all)",""),
          
          "\n\n___  PARALLEL INFO  ___\n\t",
          "Enabled:", strrep("\t",9-6), parallel.enabled,"\n\t",
          "Package used:", strrep("\t",7-4),"parallel (ver. ", packageVersion("parallel"), ")\n\t",
          "Threads utilized:", strrep("\t",5-3), parallel.threads, "\n",
          "\n"),
    file = resultOutput.info, append=FALSE)



rm(inputFileName.k1.prefix, inputFileName.k1.suffix, inputFileName.kN.prefix, inputFileName.kN.suffix, inputFileName.kN.range, inputFileName.dir.prefix)  #  Freeing up the memory by removing the variables no longer needed


runtimeStart=Sys.time()


######  Adding NT + invNT (inverse NT) occurrences if required  ######
# Check if need to double up the A/T C/G pairing occurrences
if(Add_Antisense_Occurrence==TRUE && Extended_Genome==FALSE){
  warning("Need to verify if this is correct or if it needs to be 2*occ.k1 or simply just left alone for the lone nucleotide frequency.")
  occ.k1 <- occ.k1 + occ.k1[,invNT[colnames(occ.k1)]]  #  ~Double up the total basepairs counted as the antisense bases are being added
  for(aX in 2:max(Report_K_Range)){
    curK = paste0("occ.k",aX)
    temp_occK = get(curK)   #  Assigning to another variable so get() does not need to be called multiple times
    temp_replK = temp_occK  #  Keeping the changes in another variable so multiple additions are not performed
    for(ntCol in colnames(temp_occK)){
      if(ntCol=="Other"){next}
      invSeq = paste0(invNT[strsplit(ntCol,"")[[1]]],collapse="")
      temp_replK[,ntCol] = temp_occK[,ntCol] + temp_occK[,invSeq] # Adding the inverse NTs
    }
    assign(curK, temp_replK)
  }
  rm(curK, ntCol, ntb, temp_occK, temp_replK, aX)
  
  # Ensuring that the status is marked as completed
  Add_Antisense_Occurrence = FALSE
  Extended_Genome = TRUE

  cat("\n\tDoubled the complementary basepair occurrences for obtaining the corrected abundance in the 'extended' genome (to account for potential directional bias)\n",
      file = resultOutput.info)
}


# Safety check if f*(k=1) is actually f(k=1) or not (i.e. did it properly account for the base)
# This is verified as N*(k2)/N(k1) ~ 2 [ N*(k2)/N*(k1) ~ 1] and double checked by any one of n(NT)!=n(inv(nt))
check.needs_inv_freq_accounted = matrix(0, nrow = dim(occ.k1)[1], ncol = 1, dimnames = list(rownames(occ.k1),"NULL"))
check.needs_inv_freq_accounted = cbind(check.needs_inv_freq_accounted, rowSums(occ.k1), rowSums(occ.k2))
colnames(check.needs_inv_freq_accounted) = c("need_double","N_k1","N_k2")
if(Extended_Genome==TRUE){
  for(orgName in rownames(occ.k1)){
    if(round(check.needs_inv_freq_accounted[orgName,"N_k2"]/check.needs_inv_freq_accounted[orgName,"N_k1"])!=1){
      warnMsg = paste0("'", orgName,"' will have it's f(k=1) values doubled as it is expected that the frequency was obtained in a single direction (not accounting for antisense) as it does not match the expected ratio of N_k2/N_k1")
      warning(warnMsg)
      cat(paste0(warnMsg,"\n"), file = resultOutput.info)
      check.needs_inv_freq_accounted[orgName, 1] <- 1
    }else if(DEBUG.force_complement_k1==TRUE){
      warnMsg = paste0("'", orgName,"' will have it's f(k=1) values doubled as the variable 'DEBUG.force_complement_k1' is set to TRUE!")
      warning(warnMsg)
      cat(paste0(warnMsg,"\n"), file = resultOutput.info)
      check.needs_inv_freq_accounted[orgName, 1] <- 1
    } 
  }
  rm(warnMsg)
}

asteriskFreq = ifelse(Extended_Genome==TRUE, "*", "")

######  Calculating observed frequencies  ######
cat("::: Calculating observed frequencies :::\n")

cat("\n\n___ EXPLANATIONS ___\n\n",file = resultOutput.info)


for(fi in 1:max(Report_K_Range)){
  curOK=paste0("freq.k",fi,".obs")

  ##~~   Revised method for expected freq calculations in rel. abund. section so the initial freq.k#.exp have been removed here
  #  Set up the individual base frequencies when for k=1
  if(fi==1){
    assign(curOK, occ.k1/rowSums(occ.k1))
  }else{
    occ = get(paste0("occ.k",fi)) # Obtain the variable k#
    freq.obs = occ

    for(orgName in rownames(freq.obs)){   # Loop through the rows (organism name)
      sumOcc = sum(occ[orgName,])
      for(ntSeq in colnames(freq.obs)){ # Loop through the columns (oligonucleotide groupings)

        # Observed frequency = occ(ntSeq)/sum_k(occ(nt_k))
        freq.obs[orgName, ntSeq] = freq.obs[orgName, ntSeq]/sumOcc

        }
      }
    
    # Safety check. All frequencies should add to 1
    if(round(max(rowSums(freq.obs)))!=1){ stop(paste0("Obs(f*(k=",fi,")) is NOT CORRECT when checking max freq.obs! Expected max == 1 but it is ",max(rowSums(freq.obs)))) }

    assign(curOK, freq.obs)
    
    if(output.debug.files==TRUE){
      debugOutput.freq_obs = paste0(debugOutput.freq.path, resultFile.file.prefix, "nt-freq_k",fi,".tsv")
      if(output.gzip.debug_files==TRUE){ debugOutput.freq_obs = gzfile(paste0(debugOutput.freq_obs,".gz"), "w", compression = 9, encoding = output.encoding);
      }else{debugOutput.freq_obs = file(debugOutput.freq_obs, "w", encoding = output.encoding)}
      
      if(debugOutput.freq.omit_zero_columns == TRUE){
        freq.obs = freq.obs[,intersect(
          c(sort(colnames(freq.obs)[-which(colnames(freq.obs) %in% "Other")]),"Other"),
          labels(which(colSums(freq.obs)!=0))
          )]
      }
      write.table(freq.obs, file = debugOutput.freq_obs, append=FALSE, quote=F, sep="\t", col.names = NA)
      close(debugOutput.freq_obs)
    }
    
    rm(freq.obs, occ, curOK, debugOutput.freq_obs)   #  To ensure it is not accidentally used twice in any error scenarios
  }
}
rm(fi, debugOutput.freq.path, debugOutput.freq.omit_zero_columns);


cat(paste0("\n\nThe Expected Frequency is calculated by having the product of each base of the oligomer divided by the total count/bases in k=1 (i.e. f(nt_base = oligo[i])/N(k=1)) \n",
           "If the genome is extended, then f(k=1) must be f*(k=1):\n",
           "\tf*(X) = f(X) + f(inv(X))\n",
           "This script accounts for it if it is noticed that f(k=1) was used as f*(k=2)/f(k=1) ~ 2 whereas f*(k=2)/f*(k=1) ~ 1. Theoretically, always equalizing the bases (f*(X) = f(X) + f(inv(X))) for either f/f*(k=1) should be fine as the final ratio is looked at, but for safety measure, only the suspected portion is doubled unless DEBUG.force_complement_k1 == TRUE\n\n",
           
           "Exp. Freq (oligo) = SUM_(i=1 -> k=|oligo|) [ occ(oligo[i]) / N(k=1) ]\n\t",
           "Where:\n\t\t",
           "occ = occurrence (count),\n\t\t",
           "oligo[i] = base i of oligo (going through each mononucleotide composition of oligo given),\n\t\t",
           "N(k=1) = total count/number of mononucleotide bases (word size 'k' = 1)\n\n\t\t",
           "Note 1: For 'Other', it is equal to the observed f*('Other' @ k=1)\n",
           "Note 2: For motifs containing Ns (e.g. XNZ), the sum of all the possibile frequencies is performed (e.g. f*(ANG) = sum(f*(A[ATCG]G)))\n\n"
           ),
    file = resultOutput.info)


###### Calculating Rel. Abundance and obtaining the Absolute Differences between all rows ######

registerDoParallel(pClust)  # Enable the use of parallel pthreads (%dopar%) foreach in the freq.exp portion

# Absolute relative abundance between two organisms for each kmer. Three dimensions: [x=orgName, y=orgName, z=kmer]
RelAbund.delta <- array(NA, dim = c(dim(occ.k1)[1], dim(occ.k1)[1], length(Report_K_Range)), dimnames=list(rownames(occ.k1), rownames(occ.k1), c(paste0("k",Report_K_Range))))
{
  for(kmer in Report_K_Range){
    cat(paste0("::: Calculating k=",kmer," Relative Abundance (make take several minutes, using ", parallel.threads, " pthreads) :::\n"))
    startTime=Sys.time()
    k_comb_motifs = knt_letters[1:kmer] # max(Report_K_Range)
    max_kmer_base = paste0(k_comb_motifs, collapse="")
    
    ######  Reporting the Rel. Abundance formula in info ######
    cat(paste0("Running the k=", kmer," frequency calculations\n------------  k=", kmer," (", substr(max_kmer_base, 1, kmer) ,", ", (2**kmer)-1," combinations)  -----------\n"), file = resultOutput.info)
    k_comb_motifs = derivATCG(max_kmer_base, replChar = "N", splitOddEven = FALSE)

    # Printing the combinations sorted by the number of non-N characters for the equation in info
    {
      ## Sorting in to product sizes
      # foot___ = footer [equation] notes
      motifOrder = nchar(gsub("N","",k_comb_motifs))
      footEqStart = paste0("[ f",asteriskFreq,"(")
      footEqTRUE = paste0("{  ", footEqStart)  # ...TRUE = Part of the numerator (odd/even value depends on the reporting kmer value)
      footEqFALSE=footEqTRUE                     # ...FALSE = Part of the denominator (odd/even frequencies swap as it is further derived)
      isTop=FALSE
      for(mi in max(motifOrder):1){  ##!! FOR TESTING, change range
        isTop = !isTop
        footEqPos = paste0("footEq",isTop)
        assign(footEqPos, paste0(get(footEqPos), paste0(k_comb_motifs[which(motifOrder==mi)], collapse=") f*("), collapse=""))
        assign(footEqPos, paste0(get(footEqPos), ") ]   ", footEqStart, collapse=""))
        # mi=mi-1;get(footEqPos)  ##!! TESTING DEBUG
      }

            # Removing the extra blank collapse entries
      footEqTRUE = gsub(".{7}\\($|", "",footEqTRUE); footEqFALSE = gsub(".{7}\\($", "",footEqFALSE)
      footEqTRUE = gsub("^\\s+|\\s+$", "",footEqTRUE); footEqFALSE = gsub("^\\s+|\\s+$", "",footEqFALSE)
      
      # Printing the summarized version of the formula
      baseSize = nchar(max_kmer_base)
      cat(paste0("\n",strrep("-",46),"\n---  Relative Abundance (RA) {odds ratio}  ---\n",strrep("-",46),"\n",
                  "\nRA",asteriskFreq,"(k=",baseSize,") = ",
                  switch(min(baseSize,5),"f","\U03C1","\U03B3","\U03C4","RA"),asteriskFreq,"(",max_kmer_base,") = ",
                  "{C(oligo, 2N",switch(baseSize%%2,"+1"),")|N\U03B5\U2124} \U00F7 {C(oligo, 2N",switch(1+baseSize%%2,"+1"),")|N\U03B5\U2124}",
                 "\n         = RA(f",asteriskFreq,"(k = ",ifelse(baseSize%%2,"Odd","Even"),") / f",asteriskFreq,"(k = ",ifelse(baseSize%%2,"Even","Odd"),"))",
                 "\n         = { f",asteriskFreq,"(k=",baseSize,")",ifelse(baseSize>2,paste0(" f",asteriskFreq,"(k=",baseSize-2,")"),""),ifelse(baseSize>1,paste0(" ... f",asteriskFreq,"(k=2N",switch(baseSize%%2,"+1"),")"),""),
                  "}  /  { f",asteriskFreq,"(k=",baseSize-1,")",ifelse(baseSize>3,paste0(" f",asteriskFreq,"(k=",baseSize-3,")"),""),ifelse(baseSize>1,paste0(" ... f",asteriskFreq,"(k=2N",switch(1+baseSize%%2,"+1"),")"),""),
                  " }\n"),
          file = resultOutput.info)

      # Extended equation to be printed at the end:
      cat(paste0("\n",strrep(" ",11+length(max_kmer_base)), footEqTRUE, " }\n",
                "RA",asteriskFreq,"(",max_kmer_base,") = ", strrep("-", 2+max(nchar(footEqFALSE), nchar(footEqTRUE))), "\n",
                "",strrep(" ",11+length(max_kmer_base)), footEqFALSE, " }\n\n",
                collapse=""),
          file = resultOutput.info)

      rm(footEqStart, footEqTRUE, footEqFALSE)
    }
    
    if(DEBUG>2){
        print(length(k_comb_motifs)); # Number combinations formed

        (k_comb_motifs[order(motifOrder)])
        for(mi in 1:max(motifOrder)){
          cat(paste0("\n***  dK == ",mi," --> C(", kmer, ",", mi,") = ", choose(kmer,mi),"\n"))
          print(k_comb_motifs[which(motifOrder==mi)])
      }
    }

    cat("\n")
    
    
    
    #######  CALCULATING THE REL. ABUNDANCE   #######
    
    ## DEPRECIATED:  DEBUG.use_ln_transform  [as frequency should be randing between 0~2ish, the alternative of taking large values in to account with ln-transformation should not be needed]
    
    Calculate_RelAbund_Motif <- function(orgName, iterBases, freqOBS){
      # orgName   = Name of the organism currently looped (can be the current index as well)
      # iterBases = vector of motifs used for the motif. E.g. CAT, CA, C.T, AT, C, A, T
      # freqOBS   = List of observed frequencies

      RA = 1
      for(nomDenom in 1:2){

        for(ndBase in iterBases[[nomDenom]]){
          freqKOBS = freqOBS[[nchar(ndBase)]][orgName,]
          if(grepl("[.]", ndBase)==TRUE){
            ndFreq = sum(freqKOBS[grep(ndBase, colnames(freqKOBS))])

          }else{
            ndFreq = unlist(freqKOBS[ndBase])  #  Get the f*(...) of the current base
          }

          # If in the numerator (nomDenom==1), multiply the values, otherwise divide
          RA = ifelse(nomDenom==1, RA*ndFreq, RA/ndFreq)
          names(RA) = orgName
        }
      }
      return(RA)
    }

    labelK = paste0("k",kmer)
    possMotifs = permATCG(kmer)
    RelAbund.motif <- matrix(1, nrow = dim(occ.k1)[1], ncol = 4^kmer, dimnames = list(rownames(occ.k1), possMotifs))

    listOBS=mget(paste0("freq.k", 1:kmer,".obs"))

    #$ Time consuming step, performing parallel for simultaneous organism calc. (try splitting to see if deriv parallel calc would be faster)
    #$ k=4 time (via profile) @ 47260 ~ 52210 ms
    for(kBase in possMotifs){ 
      RelAbund.motif[,kBase] = foreach(orgName = rownames(occ.k1), .combine = c, .inorder = F) %dopar% Calculate_RelAbund_Motif(orgName, derivATCG(kBase), listOBS)
    };
    
    assign(paste0("RelAbund.motif.k",kmer), RelAbund.motif)

    ## Forcing equalization among the frequencies:  forceAsterisk = (DEBUG.force_complement_k1==TRUE || check.needs_inv_freq_accounted[orgName,1]==1)
    if(DEBUG.Average_RevComp_Motifs==TRUE){
      
      Force_Avg_RevComp_Data <- function(orgName, RA_motif){
        X = RA_motif[orgName,] #Already included.
        skipMotif=c()
        for(motifBase in possMotifs){
          if(motifBase %in% skipMotif){next}
          revMotif = revcomp(motifBase)
          compMotifs = c(motifBase, revMotif)
          newVal = mean(X[compMotifs])
          X[compMotifs] = newVal
          skipMotif=c(skipMotif, revMotif)
        }
        XX = as.data.frame(t(X[order(names(X))]), row.names=orgName)
        return(XX[!(names(XX) %in% skipMotif)])
      };
      assign(paste0("RelAbund.motif.k",kmer,".orig"), RelAbund.motif)
      assign(paste0("RelAbund.motif.k",kmer), foreach(orgName = rownames(occ.k1), .combine = rbind) %dopar% Force_Avg_RevComp_Data(orgName,RelAbund.motif))
    }
    
    rm(listOBS)

    for(org1 in rownames(occ.k1)){
      for(org2 in rownames(occ.k1)){
        RelAbund.delta[org1, org2, labelK] <- sum(abs(RelAbund.motif[org1, ]-RelAbund.motif[org2, ]))
      }
    }
    RelAbund.delta[, , labelK] <- (RelAbund.delta[, , labelK]/(4**kmer)) # View(RelAbund.delta[, , labelK]*1000)
    assign(paste0("RelAbund.delta.k",kmer), RelAbund.delta[, , labelK])
    rm(org1, org2)


    #######  OUTPUT THE REL. ABUND. DIFFERENCE   #######
    stopTime = Sys.time()
    cat(paste0("\n---  k=", kmer," RUNTIME INFO  ---\n\tRun start: ", startTime, " ", format(startTime,"%Z (UTC%z)"),
               "\n\tRun stop: ", stopTime, " ", format(stopTime,"%Z (UTC%z)"),
               "\n\tRun time: ", format(stopTime-startTime),"\n"), file = resultOutput.info)
    rm(startTime, stopTime)

    # Output the delta result file
    resultOutput.Delta = paste0(resultFile.path.default, "delta_relAbund_k", kmer,".tsv")
    if(output.gzip.result_files==TRUE){ resultOutput.Delta = gzfile(paste0(resultOutput.Delta,".gz"), "w", compression = 9, encoding = output.encoding); 
    }else{resultOutput.Delta = file(resultOutput.Delta, "w", encoding = output.encoding)}
    
    write.table(RelAbund.delta[,,labelK], file = resultOutput.Delta, append=FALSE, quote=F, sep="\t", col.names = NA)
    close(resultOutput.Delta)

    if(output.debug.files==TRUE){
      debugOutput.RA_Motif = paste0(debugOutput.relAbund.motif.path, resultFile.file.prefix, "relAbund_motif_k",kmer,ifelse(DEBUG.Average_RevComp_Motifs==TRUE,"_avg-revcomp-motifs",""),".tsv")
      if(output.gzip.debug_files==TRUE){ debugOutput.RA_Motif = gzfile(paste0(debugOutput.RA_Motif,".gz"), "w", compression = 9, encoding = output.encoding); 
      }else{debugOutput.RA_Motif = file(debugOutput.RA_Motif, "w", encoding = output.encoding)}
      if(debugOutput.relAbund.motif.omit_zero_columns == TRUE){ RelAbund.motif = RelAbund.motif[,labels(which(colSums(RelAbund.motif)!=0))] }
      write.table(RelAbund.motif, file = debugOutput.RA_Motif, append=FALSE, quote=F, sep="\t", col.names = NA)
      close(debugOutput.RA_Motif)
      rm(debugOutput.RA_Motif)
    }
    rm(startTime, stopTime, resultOutput.Delta, RelAbund.motif)
  }
}; # This step may take several minutes to compute, please wait patiently
stopImplicitCluster()  #  No longer need to use cores at the moment as foreach() use is done


if(Extended_Genome==TRUE){
  cat("\n\nIt is assumed that the occurrences obtained from the genome have accounted for the complementary bases as the option has been enabled!\n",
      file = resultOutput.info)
}else if(Add_Antisense_Occurrence==TRUE){
  cat("\n\nIt is assumed that the occurrences obtained from the genome did not initially account for the complementary bases which was performed from this script as the option has been enabled!\n",
      file = resultOutput.info)
}


cat(paste0("\nf(NT) = frequency of NT = (occurrence of NT in genome X)/sum(all NT in genome X)\n",
           "f*(NT) = f(NT) + f(inv(NT)) = frequency of the 'extended' genome (both sense and antisense present). Inv = inverse (complementary NT base)"
           ),
    file = resultOutput.info)




####  PCA - Principal Component Analysis  ####
# For more info, see http://sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp
#                    https://www.sartorius.com/en/knowledge/science-snippets/what-is-principal-component-analysis-pca-and-how-it-is-used-507186
#                    https://ml-explained.com/blog/principal-component-analysis-explained

if(options.Output_PCA==TRUE){
  cat("::: Running PCA analysis :::\n")
  startTime=Sys.time()

  # Preparing PCA output type
  PCA.Image_Output = tolower(PCA.Image_Output)
  if(library("grDevices", character.only=TRUE, logical.return=TRUE)==FALSE){  # Should already be installed by default but just in case... Version tested: 4.1.0
    install.packages(grDevices) # Version tested: 2.0.0
    require(grDevices)
  }
  
  checkPath = dirname(PCA.Output.Prefix); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)

  ### IMPORTANT FUNCTIONS ###
  plotPCA <- function(Data.All,                          # Data = The data to be analyzed (e.g. data.raw[,3:19]) 
                      assignPCA2var = "",                # A string with the variable to save the prcomp() data is saved
                      Title         = "",                # Title = The Title part of the dataset analyzed for labeling purposes (e.g. "IMP+IRO+IE")
                      SetName       = dataset.id,        # Dataset name = Part of the filename to help differentiate
                      outputImage   = PCA.Image_Output,  # makePNG = If TRUE, make Image. If FALSE, plot directly
                      repelText     = TRUE,              # repelText = For the PCA graph, repel text or not
                      addEllip      = PCA.Add_Ellipsis,  # addEllip = add ellipses to the PCA plot
                      ellipFactor   = "",                # ellipFactor = The factor for the groupings to be considered for the ellipsces (only if addEllip is set to TRUE)
                      ellipLevel    = PCA.Ellipsis.CI,   # Confidence interval
                      ellipTitle    = "Class",           # Title for the symbol legend
                      plotScree     = PCA.Scree_Plot,    # plotScree = Output Scree (eigenvalue) plot as well
                      plotOnlyScree = FALSE,             # If plotScree is TRUE, then use this to stop displaying the regular PCA directly after
                      Image.width   = PCA.Width,
                      Image.height  = PCA.Height,
                      centerData    = TRUE,
                      scaleData     = FALSE
  ){
    if(!exists(PCA.Function) || (PCA.Function!="prcomp" && PCA.Function!="princomp")){ stop("PCA.Function not correctly defined (use either prcomp [recommended] or princomp!)") }
    Data.All.pca<- get(PCA.Function)(Data.All, center = centerData, scale = scaleData)
    if(assignPCA2var!=""){assign(assignPCA2var, Data.All.pca, envir=.GlobalEnv)}
    if(outputImage=="pdf"||outputImage=="svg"){
      Image.width = round(Image.width/Pixels_Per_Inch)
      Image.height = round(Image.height/Pixels_Per_Inch)
    }
    if(plotScree==TRUE){
      if(outputImage!=""){
        if(!exists(outputImage)){ stop2("ERROR: image function does not exist (see PCA.Image_Output)") }
        get(ifelse(outputImage=="svg","svglite",outputImage))(gsub(Regex_Remove_Filename_Chars, "", 
                              gsub(" ", "_", paste0(PCA.Output.Prefix, "_PCA-Scree-plot_", SetName, "_", Title,".",outputImage)), 
                          ignore.case = TRUE), width=Image.width, height=Image.height)
      }
      print(fviz_eig(Data.All.pca, title = paste(Title, " (",SetName, ") SCREE PLOT")))
      if(outputImage!=""){ dev.off() }
      if(plotOnlyScree==TRUE){return(NULL)}
    }

    if(outputImage!=""){
      if(!exists(outputImage)){ stop2("ERROR: image function does not exist (see PCA.Image_Output)") }
      get(ifelse(outputImage=="svg","svglite",outputImage))(gsub(Regex_Remove_Filename_Chars, "", 
                            gsub(" ", "_", paste0(PCA.Output.Prefix, "_PCA-plot_", SetName, "_", Title,".",outputImage)), 
                        ignore.case = TRUE), width=Image.width, height=Image.height)
      # , format(Sys.time(),"_%y%m%d_%H%M%S")
    }
    if(addEllip==FALSE){
      if(length(ellipFactor)<=1){
        # Grouping NOT passed along
        print(fviz_pca_var(
          Data.All.pca,
          title = paste0(Title, " (",SetName, ")"),
          geom = c("text","point"),
          
          ##habillage = ellipFactor, # Group by grouping factor given
          ##show.legend=F, # Gets rid of the text that appears in the legend 
          ##legend.title = ellipTitle,
          col.var = "contrib", # Color by contributions to the PC
          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
          repel = repelText     # Avoid text overlapping
        ))
      }else{
        # Grouping passed along
        print(fviz_pca_var(
          Data.All.pca,
          title = paste0(Title, " (",SetName, ")"),
          geom = c("text","point"),
          
          habillage = ellipFactor, # Group by grouping factor given
          show.legend=F, # Gets rid of the text that appears in the legend 
          legend.title = ellipTitle,
          ## col.var = "contrib", # Color by contributions to the PC
          ## gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
          repel = repelText     # Avoid text overlapping
        ))
      }
    }
    else{
      ##groups <- as.factor(ellipFactor)
      print(fviz_pca_var(Data.All.pca,
                         title = paste0(Title, " (",SetName, ")"),
                         geom = c("text","point"),
                         repel = repelText,     # Avoid text overlapping
                         ##col.var = groups, # Color by contributions to the PC
                         habillage=ellipFactor,
                         addEllipses = TRUE, # Concentration ellipses
                         ellipse.type = "confidence",
                         ellipse.level=ellipLevel,
                         show.legend=F, # Gets rid of the text that appears in the legend 
                         legend.title = ellipTitle
      ))
    }
    
    if(outputImage!=""){dev.off()}
  }
  
  # === EXAMPLE OF PLOTTING THE Title (VARIABLE data.all = COLUMN data.all DISPLAYED) ===
  for(kmer in Report_K_Range){
    freq = t(get(paste0("RelAbund.motif.k",kmer))) # rownames(freq) = colnames(occ.k1)[1];
  
    # print(fviz_pca_var(prcomp(t(ref.freq[,-1]),center=T,scale=F),geom = c("point","text"),repel=T,habillage=ref.freq[,1]))
    plotPCA(freq,
            Title = paste0("k=", kmer, " motif PCA"), 
            assignPCA2var="PCA.Plot_Data",
            plotScree = TRUE # Prints both Scree (eigenvalue) and regular PCA component 1 vs 2 plots
            ) 
    # If you get an error stating Rpp doesn't have 'Rcpp_precious_remove', reinstall Rcpp with either:
    #   update.packages()
    #   or:  install.packages("Rcpp") # Must be version >= 1.0.7
    
    
    # Summary
    outFile.tmpName = gsub("PCA/",paste0("PCA/k",kmer,"/summary/"), PCA.Output.Prefix); checkPath = dirname(outFile.tmpName); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
    outFile.temp = gzfile(gsub(Regex_Remove_Filename_Chars, "", 
                               paste0(outFile.tmpName,"PCA-data_summary_", dataset.id, "_k", kmer, ".tsv.gz")),
                          "w", compression = 9, encoding = output.encoding) 
    write.table(summary(PCA.Plot_Data)$importance, outFile.temp, append=FALSE, quote=F, sep="\t", col.names = NA)
    close(outFile.temp)

    # Eigenvalue (relates to proportion)    
    outFile.tmpName = gsub("PCA/",paste0("PCA/k",kmer,"/eigenvalue/"), PCA.Output.Prefix); checkPath = dirname(outFile.tmpName); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
    outFile.temp = gzfile(gsub(Regex_Remove_Filename_Chars, "", 
                               paste0(outFile.tmpName,"PCA-data_eigenvalue_", dataset.id, "_k", kmer, ".tsv.gz")),
                          "w", compression = 9, encoding = output.encoding) 
    write.table(get_eigenvalue(PCA.Plot_Data), outFile.temp, append=FALSE, quote=F, sep="\t", col.names = NA)
    close(outFile.temp)
    
    # Contribution of the motifs (ind)
    outFile.tmpName = gsub("PCA/",paste0("PCA/k",kmer,"/ind_contrib/"), PCA.Output.Prefix); checkPath = dirname(outFile.tmpName); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
    outFile.temp = gzfile(gsub(Regex_Remove_Filename_Chars, "", 
                               paste0(outFile.tmpName,"PCA-data_motif-contribution_", dataset.id, "_k", kmer, ".tsv.gz")),
                          "w", compression = 9, encoding = output.encoding) 
    write.table(get_pca_ind(PCA.Plot_Data)$contrib, outFile.temp, append=FALSE, quote=F, sep="\t", col.names = NA)
    close(outFile.temp)
    
    # Contribution of the individuals (var)
    ## get_pca_var(PCA.Plot_Data) <-- species contrib.
    ## get_pca_ind(PCA.Plot_Data) <-- motif contributions
    outFile.tmpName = gsub("PCA/",paste0("PCA/k",kmer,"/var_contrib/"), PCA.Output.Prefix); checkPath = dirname(outFile.tmpName); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
    outFile.temp = gzfile(gsub(Regex_Remove_Filename_Chars, "", 
                               paste0(outFile.tmpName,"PCA-data_motif-contribution_", dataset.id, "_k", kmer, ".tsv.gz")),
                          "w", compression = 9, encoding = output.encoding) 
    write.table(get_pca_var(PCA.Plot_Data)$contrib, outFile.temp, append=FALSE, quote=F, sep="\t", col.names = NA)
    close(outFile.temp)
    
    
    rm(PCA.Plot_Data, outFile.tmpName)

  }
  ## rm(PCA.Add_Ellipsis, PCA.Ellipsis.CI, PCA.Height, PCA.Image_Output, PCA.Output.Prefix, PCA.Plot_Data, PCA.Scree_Plot, PCA.Width)
  
  stopTime = Sys.time()
  cat(paste0("\n\n___  PCA RUNTIME INFO  ___\n\tRun start: ", startTime, " ", format(startTime,"%Z (UTC%z)"),
             "\n\tRun stop: ", stopTime, " ", format(stopTime,"%Z (UTC%z)"),
             "\n\tRun time: ", format(stopTime-startTime),"\n"), file = resultOutput.info)
  rm(startTime, stopTime)
  
} # PCA section


####  NJ TREE ####

# Phylogenetic view
# APE version 5.1 (A tutorial is available at http://r-eco-evo.blogspot.com/2007/09/neighbor-joining-tree-with-ape.html)

#  WHAT IS NEEDED
#  ---
#  A square matrix file needs to be provided with the values you want to plot out.
#     (There is a small conversion part in this file if you need to transform your similiarity matrix in to a
#      distance/dissimilarity matrix or vice-versa if needed)
#     Ensure that the appropriate entries are spaced out correctly

# Checks if APE is installed. If not, install it.
# Load the k=2 table


if(options.Output_NJ_Tree==TRUE){
  cat("::: Plotting Neighbour-Joining Tree using the Absolute Difference (Delta) Relative Abundance values between the different organisms :::\n")
  startTime = Sys.time()

  # Preparing Tree image output (PHYLIP bootstrap is done separately)
  if(Tree.Image_Output!=""){
    if(Tree.Image_Output=="pdf"||Tree.Image_Output=="svg"){
      PCA.Width  = round(PCA.Width/Pixels_Per_Inch)  # Since PCA is done, we are reusing the same variable for the Tree (if defined to output an image)
      PCA.Height = round(PCA.Height/Pixels_Per_Inch)
    }
    Tree.Image_Output = tolower(Tree.Image_Output)
    if(library("grDevices", character.only=TRUE, logical.return=TRUE)==FALSE){  # Should already be installed by default but just in case... Version tested: 4.1.0
      install.packages(grDevices) # Version tested: 2.0.0
      require(grDevices)
    }
    get(ifelse(Tree.Image_Output=="svg", "svglite", Tree.Image_Output))(
      gsub(Regex_Remove_Filename_Chars, "", 
           gsub(" ", "_", paste0(Tree.Output.Prefix, "motif_heatmap_k", kmer,".",HeatMotif.Image_Output)), 
           ignore.case = TRUE), width=PCA.Width, height=PCA.Height)
  }
  
  
  checkPath = dirname(Tree.Output.Prefix); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)


  njTree <- function(
    Data,           # The data set to analyze
    rootTo      = Tree.Root_To,      # Row name of what to root to (leave blank for unrooted)
    Bootstrap   = Tree.Bootstrap, # The number of bootstraps needed
    saveTree=""     # The name of the file to be saved to. Appends _unrooted.phb/_rooted-[root].phb
  ){
    njtree.data <- Data; head(njtree.data);dim(njtree.data)  #  View the first 5 data lines to ensure everything is in correctly
    
    # Producing a neighbor-joining tree with a thousand bootstrap support
    #    To read more about the nj() function from the APE package, run: ?nj
    #    A function was defined as it needs to be used in the bootstrap function below
    
    njtree.func <- function(x) nj(dist(x))
    # Alternatively, you could use hclust() for the clustering instead of 
    # as.dist() which is used for the computation of the distance matrix
    
    njtree.phylo.unroot <- njtree.func(njtree.data)
    njtree.isRooted <- FALSE 
    if(rootTo==""){
      njtree.phylo <- njtree.phylo.unroot
    }else{
      njtree.phylo.root <- root(phy = njtree.phylo.unroot, outgroup = rootTo)
      njtree.phylo <- njtree.phylo.root
      njtree.isRooted <- TRUE
    }
    
    # Set which tree to use (here, we went with the rooted tree)
    
    # Producing a 1,000 bootstrap phylogeny.
    #   phy  = the phylogram object which was generated above
    #   x    = the data used for generating the phylogram
    #   FUN  = the function used to estimate the phylogram object (phy)
    #   B    = the number of bootstraps to be performed
    # rooted = If the tree should be rooted or not
    # trees  = TRUE to return the tree/phylo object (if set to false, it returns just the bootstrap values,
    #          which is nice to see if the bootstraps are consistent. See ?boot.phylo for the example)
    
    # Run the following commented lines to see if the bootstrap support changes:
    # for(i in 1:5){print(boot.phylo(njtree.phylo, njtree.data, njtree.func, B=1000, quiet=TRUE))}
    
    njtree.boot <- boot.phylo(phy = njtree.phylo, x = njtree.data, FUN = njtree.func,
                              B=Bootstrap, trees=TRUE, rooted=njtree.isRooted)
    
    njtree.clades<-prop.clades(njtree.phylo, njtree.boot$trees, rooted=njtree.isRooted)
    
    njtree.phylo$node.label <- njtree.clades   #  Place the bootstrap values on the nodes. Obtained from https://www.biostars.org/p/9511/
    
    # Display  the rooted bootstrap NJ tree
    plot(njtree.phylo, show.node.label = TRUE)
    if(Tree.Image_Output!=""){dev.off()}
    
    # nodelabels(njtree.clades)
    # drawSupportOnEdges(njtree.boot$BP) # If you wanted the bootstrap values on the edges rather than the nodes
    
    
    # Export the bootstrap NJ tree
    njtree.phylo$node.label[is.na(njtree.phylo$node.label)]<-""   # Convert NA node labels to blanks
    
    # Exports as a Newick tree. For a NEXUS tree format, please see ?write.nexus BUT note that NEXUS format has issues exporting the node labels (as mentioned in https://www.researchgate.net/post/SOLVED_How_do_you_export_bootstrap_node_support_in_Rs_ape_package)
    if(saveTree!=""){
      njtree.filename = gsub(Regex_Remove_Filename_Chars, "", 
                             paste0(
                               saveTree, # "njtree_",
                               ifelse(njtree.isRooted,paste0("_rooted-",rootTo),"_unrooted"), 
                               ".phb"
                             ), 
                             ignore.case = TRUE)
      # print(njtree.filename)
      write.tree(
        phy = njtree.phylo,
        file = njtree.filename
      )
    }
  }
  
  #Tree.Output.Prefix    #    Appends _unrooted.phb/_rooted-[root].phb depending on the Tree.Root_To value
  ## par(mfrow=c(1,length(Report_K_Range)))  # Show multiple images at once (row, col)
  for(kmer in Report_K_Range){
    njTree(get(paste0("RelAbund.delta.k",kmer)), 
           saveTree=paste0(Tree.Output.Prefix,Tree.Bootstrap,"-boot_NJ-tree_Rel.Abund_k",kmer) # saveTree suffix with will be appended with '_unrooted.phb/_rooted-[root].phb' in the function
           )
    
    title(paste0("Absolute Difference Relative Abundance (k=", kmer, ")"))  # Generated average delta values
  }
  ## par(mfrow=c(1,1)) # Set back to showing just one image output
  
  stopTime = Sys.time()
  cat(paste0("\n\n___  NJ TREE RUNTIME INFO  ___\n\tRun start: ", startTime, " ", format(startTime,"%Z (UTC%z)"),
             "\n\tRun stop: ", stopTime, " ", format(stopTime,"%Z (UTC%z)"),
             "\n\tRun time: ", format(stopTime-startTime),"\n"), file = resultOutput.info)
  rm(startTime, stopTime)
  
} # NJ tree section




####  FREQUENCY HEATMAP (MOTIF)  ####
# Graphing options: https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/

##### Rel. Abund. Motif Heatmap ####

if(options.Output_Heatmap_Motifs==TRUE){
  cat("::: Generating Heatmap (Motif) :::\n")
  startTime = Sys.time()
  ## Heatmap of the doublet frequency
  if(HeatMotif.Grouping_File != ""){
    tmp = read.delim(HeatMotif.Grouping_File); #View(tmp)
    HeatMotif.Group_Name   <- tmp[,-1]; rownames(HeatMotif.Group_Name)<-tmp[,1]  #DEBUG: gsub("^([^_]+).+$","\\1",rownames(occ.k1)) <-- Use the genus as the grouping factor
    #HeatMotif.Group_Factor <- foreach(x = colnames(HeatMotif.Group_Name), .inorder = FALSE, .combine = list) %do% { unique(as.factor(HeatMotif.Group_Name[,x]))} #apply(HeatMotif.Group_Name, 2, )
    HeatMotif.Group_Colour <- foreach(x = colnames(HeatMotif.Group_Name), .inorder = FALSE, .combine = c) %do% {
                                  u <- unique(HeatMotif.Group_Name[,x])
                                  u = u[order(u)]
                                  n <- HeatMotif.Group_Name[,x]
                                  r=rainbow(n=length(u))#[match(n,u)]
                                  names(r) = u#rownames(HeatMotif.Group_Name)
                                  l=list(r)
                                  names(l)=x
                                  l
                                  }
    # HeatMotif.Group_Colour <- foreach(x = colnames(HeatMotif.Group_Name), .inorder = FALSE, .combine = cbind) %do% {
    #   u <- unique(HeatMotif.Group_Name[,x])
    #   n <- HeatMotif.Group_Name[,x]
    #   r=rainbow(n=length(u))[match(n,u)]
    #   names(r) = rownames(HeatMotif.Group_Name)
    #   r
    # } #apply(HeatMotif.Group_Name, 2, )

  }
  
  # Show just the k=2 rel. abund. ratio amongst the isolates (as the others will add a lot more data)
  checkPath = dirname(HeatMotif.Output.Prefix); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
  
  
  for(kmer in Report_K_Range){ 
    
    # Generate the motif combinations for the labels
    HeatMotif.Data <- get(paste0("RelAbund.motif.k",kmer))
    ignore.RevComp <- c() # Just a temp. list to skip the revcomp to reduce the combinations

    if(HeatMotif.Avg_RevComp_Motifs == TRUE){
      for(curMotif in permATCG(kmer)){
        if(curMotif %in% ignore.RevComp){  #   If RevComp already present, continue on.
          next
        }
        revMotif = revcomp(curMotif)
        ignore.RevComp = c(ignore.RevComp, revMotif)
        if(revMotif != curMotif){
          ignore.NewVector = c(curMotif, revcomp(curMotif))
          HeatMotif.Data = cbind(
            HeatMotif.Data[,-which(colnames(HeatMotif.Data) %in% ignore.NewVector)],
            rowMeans(HeatMotif.Data[,ignore.NewVector])
          )
          colnames(HeatMotif.Data)[dim(HeatMotif.Data)[2]] = paste0(curMotif, ":", revcomp(curMotif))  # Concatenate the current motif to that name and continue on
        }
      }
    } # Combining reverse complemented motifs
    ### Do not need to reorder by column name, but if you want to: HeatMotif.Data = HeatMotif.Data[,order(colnames(HeatMotif.Data))] # Order by motif
    rm(ignore.RevComp)

  # Adjusting the colouring such that white is in the middle if possible (and kmer=2 colour range is set up specifically for the distance rating by Karlin et al.)
  { #if(!exists("rampCols")){
    ignore.breaks <- seq(from=min(range(HeatMotif.Data)), to=max(range(HeatMotif.Data)), length.out=100)
    
    if(kmer==2){
      # Colour range specifically for the genomic signature data
      # Genomic signature breaks:  see legend of Table on Pg 191 of Karlin et al. 1998 or in text of Karlin and Mrazek 1997, the latter which states "We distinguish extremes of dinucleotide relative abundances as follows: extremely high, symbolically 111, r*XY $ 1.50; very high, 11, 1.30 # r*XY , 1.50; significantly high, 1, 1.23 # r*XY , 1.30; marginally high, (1), 1.20 # r*XY , 1.23; extremely low, 222, r*XY # 0.50; very low, 22, 0.50 , r*XY # 0.70; significantly low, 2, 0.70 , r*XY # 0.78; marginally low, (2), 0.78 , r*XY # 0.81."
      #                           1   2     3     4    5  |  6     7     8     9   10
      HeatMotif.Legend_Breaks = c(0,0.50, 0.70, 0.78, 0.81, 1.20, 1.23, 1.30, 1.50,2)
      
      ignore.breakpoint <- c(
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[2])), #- 0.50)), # 1 rho*(XY)<=0.50  --- extremely low
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[3])), #- 0.70)), # 2  (0.50, 0.70]   --  very low
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[4])), #- 0.78)), # 3  (0.70, 0.78]    -  significantly low
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[5])), #- 0.81)), # 4- (0.78, 0.81]   (-) marginally low
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[6])), #- 1.20)), # -5 [1.20, 1.23)   (+) marginally high
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[7])), #- 1.23)), # 6  [1.23, 1.30)    +  significantly high
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[8])), #- 1.30)), # 7  [1.30, 1.50)   ++  very high
        which.min(abs(ignore.breaks - HeatMotif.Legend_Breaks[9])) #- 1.50))   # 8 rho*(XY)>=1.50  +++ extremely high
      )
      rampCols <- c(colorRampPalette(c("darkblue", "blue"))(ignore.breakpoint[1]),
                  colorRampPalette(c("blue", "#90D2EC"))(ignore.breakpoint[2]-ignore.breakpoint[1]),
                  colorRampPalette(c("#90D2EC", "#ECFFFB"))(ignore.breakpoint[3]-ignore.breakpoint[2]),
                  colorRampPalette(c("#ECFFFB", "white", "#FFE7BE"))(ignore.breakpoint[5]-ignore.breakpoint[3]), ##rep.int("#FFFFFF", ignore.breakpoint[4]-ignore.breakpoint[3]),
                  colorRampPalette(c("#FFE7BE", "pink"))(ignore.breakpoint[6]-ignore.breakpoint[5]),
                  colorRampPalette(c("pink", "red"))(ignore.breakpoint[7]-ignore.breakpoint[6]),
                  colorRampPalette(c("red", "darkred"))(100-ignore.breakpoint[8]+1)
                  )
    }else{
      HeatMotif.Legend_Breaks = seq(0,2,0.25)
      ignore.breakpoint <- c(
        which.min(abs(ignore.breaks - 0.99)),
        which.min(abs(ignore.breaks - 1.01))
      )
      rampCols <- c(colorRampPalette(c("darkblue", "blue", "white"))(ignore.breakpoint[1]),
                    #colorRampPalette(c("#FFFF00", "#FFFFFF", "#00FFFF"))(ignore.breakpoint[2]-ignore.breakpoint[1]),
                    rep.int("#FFFFFF", ignore.breakpoint[2]-ignore.breakpoint[1]),
                    colorRampPalette(c("white", "pink", "red", "darkred"))(100-ignore.breakpoint[2]+1)#(100-(midpoint+1))
      )
      
    }
  }
  #ignore.Group_Order <- order(HeatMotif.Group_Factor) #as.factor(classRef[ogName,"family"])); # Reorder by family first
  HeatMotif.Image_Output = tolower(HeatMotif.Image_Output)
  tmp.W = HeatMotif.Width; tmp.H = HeatMotif.Height
  if(HeatMotif.Image_Output!=""){
    if(!exists(HeatMotif.Image_Output)){ stop2("ERROR: image function does not exist (see HeatMotif.Image_Output)") }
    if(HeatMotif.Image_Output=="svg" || HeatMotif.Image_Output=="pdf"){  # Adjusting values from pixels to inches (roughly)
      HeatMotif.Width =  round(tmp.W/Pixels_Per_Inch) # ~14"
      HeatMotif.Height = round(tmp.H/Pixels_Per_Inch) # ~9"
      # Using svglite() for more appropriate standards than the default svg()
    }
    get(ifelse(HeatMotif.Image_Output=="svg", "svglite", HeatMotif.Image_Output))(gsub(Regex_Remove_Filename_Chars, "", 
                          gsub(" ", "_", paste0(HeatMotif.Output.Prefix, "motif_heatmap_k", kmer,".",HeatMotif.Image_Output)), 
                          ignore.case = TRUE), width=HeatMotif.Width, height=HeatMotif.Height)
  }

  
  # For group colouring, followed Zhigang Lu's blog post: https://zhiganglu.com/post/pheatmap_change_annotation_colors/
  dimH <- dim(HeatMotif.Data)[1]
  dimW <- dim(HeatMotif.Data)[2]
  X <- pheatmap(HeatMotif.Data, col = rampCols, na_col="#CCCCCC",
           scale="none", breaks=ignore.breaks,
           legend_breaks = HeatMotif.Legend_Breaks, # Genomic signature breaks:  see legend of Table on Pg 191 of Karlin et al. 1998 or in text of Karlin and Mrazek 1997, the latter which states "We distinguish extremes of dinucleotide relative abundances as follows: extremely high, symbolically 111, r*XY $ 1.50; very high, 11, 1.30 # r*XY , 1.50; significantly high, 1, 1.23 # r*XY , 1.30; marginally high, (1), 1.20 # r*XY , 1.23; extremely low, 222, r*XY # 0.50; very low, 22, 0.50 , r*XY # 0.70; significantly low, 2, 0.70 , r*XY # 0.78; marginally low, (2), 0.78 , r*XY # 0.81."
           cluster_rows = TRUE, annotation_row = HeatMotif.Group_Name, annotation_colors = HeatMotif.Group_Colour, annotation_names_row=T,
           cluster_cols = TRUE,
           cellwidth=round(max(500,HeatMotif.Width-450-(25*kmer))/dimW),     # Default min: 757
           cellheight=round(max(500,HeatMotif.Height-200)/dimH), # Default min: 539
           display_numbers = (kmer<4), number_color = "black", fontsize_number = ifelse(kmer==3,round(dimW/5),dimW), #k3 ~ 32/5 = 6 
           fontsize_col = (12-round(kmer*2)), main = paste0("k=",kmer, " relative abundance heatmap"), 
           angle_col = "90"
           )
   HeatMotif.Width = tmp.W; HeatMotif.Height = tmp.H
  if(HeatMotif.Image_Output!=""){ dev.off() }

  rm(ignore.breaks, ignore.breakpoint, tmp.H, tmp.W)
  
  }
  if(HeatMotif.Image_Output!=""){ dev.off() } # Just in case, make sure it is turned off
  
  stopTime = Sys.time()
  cat(paste0("\n\n___  HEATMAP (MOTIF) RUNTIME INFO  ___\n\tRun start: ", startTime, " ", format(startTime,"%Z (UTC%z)"),
             "\n\tRun stop: ", stopTime, " ", format(stopTime,"%Z (UTC%z)"),
             "\n\tRun time: ", format(stopTime-startTime),"\n"), file = resultOutput.info)
  rm(startTime, stopTime)
} # Motif heatmap section

####  FREQUENCY HEATMAP (DELTA REL. ABUND.)  ####
# Graphing options: https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/

##### Rel. Abund. Delta Heatmap (Absol. Diff. between two organisms) ####
if(options.Output_Heatmap_Delta==TRUE){
  cat("::: Generating Heatmap (Delta/Absolute diff. rel. abund.) :::\n")
  startTime = Sys.time()
  
  ## Heatmap of the doublet frequency
  if(HeatDelta.Grouping_File != ""){
    tmp = read.delim(HeatDelta.Grouping_File); #View(tmp)
    HeatDelta.Group_Name   <- tmp[,-1]; rownames(HeatDelta.Group_Name)<-tmp[,1]  #DEBUG: gsub("^([^_]+).+$","\\1",rownames(occ.k1)) <-- Use the genus as the grouping factor
    #HeatDelta.Group_Factor <- foreach(x = colnames(HeatDelta.Group_Name), .inorder = FALSE, .combine = list) %do% { unique(as.factor(HeatDelta.Group_Name[,x]))} #apply(HeatDelta.Group_Name, 2, )
    HeatDelta.Group_Colour <- foreach(x = colnames(HeatDelta.Group_Name), .inorder = FALSE, .combine = c) %do% {
      u <- unique(HeatDelta.Group_Name[,x])
      u = u[order(u)]
      n <- HeatDelta.Group_Name[,x]
      r=rainbow(n=length(u))#[match(n,u)]
      names(r) = u#rownames(HeatDelta.Group_Name)
      l=list(r)
      names(l)=x
      l
    }
    # HeatDelta.Group_Colour <- foreach(x = colnames(HeatDelta.Group_Name), .inorder = FALSE, .combine = cbind) %do% {
    #   u <- unique(HeatDelta.Group_Name[,x])
    #   n <- HeatDelta.Group_Name[,x]
    #   r=rainbow(n=length(u))[match(n,u)]
    #   names(r) = rownames(HeatDelta.Group_Name)
    #   r
    # } #apply(HeatDelta.Group_Name, 2, )
    
    # HeatDelta.Group_Colour <- gsub("FFBF00","D67702",HeatDelta.Group_Colour)  # Replace all the bright orange with a darker shade
  }
  
  # Show just the k=2 rel. abund. ratio amongst the isolates (as the others will add a lot more data)
  checkPath = dirname(HeatDelta.Output.Prefix); if(!dir.exists(checkPath)){ dir.create(checkPath, recursive = TRUE) }; rm(checkPath)
  
  
  for(kmer in Report_K_Range){ 
    
    HeatDelta.Data <- get(paste0("RelAbund.delta.k",kmer)); # View(HeatDelta.Data)
    
    if((HeatDelta.Cluster_Data == FALSE) && (length(HeatDelta.Order_Data_Col)>1 || HeatDelta.Order_Data_Col >= 0) && (HeatDelta.Grouping_File != "")){
      if((length(HeatDelta.Order_Data_Col)==1) && (HeatDelta.Order_Data_Col==0)){ HeatDelta.Order_Data_Col = dim(HeatDelta.Group_Name)[2]:1 }
      sortedDelta.order = do.call(what=order,args=HeatDelta.Group_Name[,dim(HeatDelta.Group_Name)[2]:1])
      sortedDelta.names = rownames(HeatDelta.Group_Name)[sortedDelta.order]; if(DEBUG>3){head(HeatDelta.Group_Name[sortedDelta.names,])}
      HeatDelta.Data = HeatDelta.Data[sortedDelta.names,sortedDelta.names]; if(DEBUG>3){head(HeatDelta.Group_Name[rownames(HeatDelta.Data),])}
      rm(sortedDelta.order)
    }

    # Adjusting the colouring such that white is in the middle if possible (and kmer=2 colour range is set up specifically for the distance rating by Karlin et al.)
    { 
      ignore.breaks <- seq(from=min(range(HeatDelta.Data)), to=max(range(HeatDelta.Data)), length.out=100)
      
      if(kmer==2){
        # Colour range specifically for the genomic signature data
        ### Genomic signature differences (delta) breaks:
        ## EUKARYOTIC BREAKS: :  see in text of Karlin and Mrazek 1997: extremely high, symbolically 111, r*XY $ 1.50; very high, 11, 1.30 # r*XY , 1.50; significantly high, 1, 1.23 # r*XY , 1.30; marginally high, (1), 1.20 # r*XY , 1.23; extremely low, 222, r*XY # 0.50; very low, 22, 0.50 , r*XY # 0.70; significantly low, 2, 0.70 , r*XY # 0.78; marginally low, (2), 0.78 , r*XY # 0.81."
        ##                          1   2     3     4    5     6     7      8    9   #    1     2      3     4    5     6     7     8     9   10    11     12    13
        HeatDelta.Legend_Breaks = c(0,0.018,0.030,0.050,0.075,0.115,0.150,0.200,0.250) #c(0.018,0.020,0.030,0.035,0.050,0.55,0.075,0.080,0.115,0.120,0.150,0.160,0.200)
        HeatDelta.Legend_Labels = c("random/extremely close","very close","close","moderately close", "weakly close", "weakly distant", "distant", "very distant", "")

        ignore.breakpoint <- c(
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[2])), #- 0.018,  # 1 delta*(G1,G2)<=0.018  random
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[3])), #- 0.030,  # 2  [0.020, 0.030]  very close
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[4])), #- 0.050,  # 3  [0.035, 0.050]  close
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[5])), #- 0.075,  # 4  [0.055, 0.075]  moderately similar
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[6])), #- 0.115,  # 5  [0.080, 0.115]  weakly similar
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[7])), #- 0.150,  # 6  [0.120, 0.150]  distantly similar
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[8])), #- 0.200,  # 7  [0.160, 0.200)  distant
          which.min(abs(ignore.breaks - HeatDelta.Legend_Breaks[9]))  #- 0.2[5]0 # 8 delta*(XY)>=0.200  very distant
        )

        # Colour list: https://r-charts.com/colors/
        rampCols <- c(colorRampPalette(c("darkred", "red"))(ignore.breakpoint[1]),                 # Random (i.e. same species)
                      colorRampPalette(c("red","lightcoral"))(ignore.breakpoint[2]-ignore.breakpoint[1]), # very close
                      colorRampPalette(c("lightcoral", "pink"))(ignore.breakpoint[3]-ignore.breakpoint[2]),      # close
                      colorRampPalette(c("pink", "lavenderblush"))(ignore.breakpoint[4]-ignore.breakpoint[3]),  # moderately similar
                      colorRampPalette(c("lavenderblush","white"))(ignore.breakpoint[5]-ignore.breakpoint[4]),# *weakly similar
                      colorRampPalette(c("white", "lightblue"))(ignore.breakpoint[6]-ignore.breakpoint[5]),# *distantly similar
                      colorRampPalette(c("lightblue", "blue"))(ignore.breakpoint[7]-ignore.breakpoint[6]),  # distant
                      colorRampPalette(c("blue", "darkblue"))(ignore.breakpoint[8]-ignore.breakpoint[7]), # very distant
                      colorRampPalette(c("darkblue", "black"))(100-ignore.breakpoint[8]) # (very distant)
        )
      }else{
        HeatDelta.Legend_Breaks = seq(0,round(0.013+max(ignore.breaks),4),0.025)
        ignore.breakpoint = HeatDelta.Legend_Breaks
        HeatDelta.Legend_Labels = HeatDelta.Legend_Breaks

        rampCols <- colorRampPalette(c("darkred", "red", "pink","lightblue","blue","darkblue"))(100)

      }
      
    }

    #ignore.Group_Order <- order(HeatDelta.Group_Factor) #as.factor(classRef[ogName,"family"])); # Reorder by family first
    HeatDelta.Image_Output = tolower(HeatDelta.Image_Output)
    tmp.W = HeatDelta.Width; tmp.H = HeatDelta.Height
    if(HeatDelta.Image_Output!=""){
      if(!exists(HeatDelta.Image_Output)){ stop2("ERROR: image function does not exist (see HeatDelta.Image_Output)") }
      if(HeatDelta.Image_Output=="svg" || HeatDelta.Image_Output=="pdf"){  # Adjusting values from pixels to inches (roughly)
        HeatDelta.Width =  round(tmp.W/Pixels_Per_Inch) # ~14"
        HeatDelta.Height = round(tmp.H/Pixels_Per_Inch) # ~9"
        # Using svglite() for more appropriate standards than the default svg()
      }
      get(ifelse(HeatDelta.Image_Output=="svg", "svglite", HeatDelta.Image_Output))(
        gsub(Regex_Remove_Filename_Chars, "", 
             gsub(" ", "_", paste0(HeatDelta.Output.Prefix, "delta_heatmap_",ifelse(HeatDelta.Cluster_Data==TRUE,"clustered_",""),"k", kmer,".",HeatDelta.Image_Output)), 
             ignore.case = TRUE), width=HeatDelta.Width, height=HeatDelta.Height)
    }
    
    
    # For group colouring, followed Zhigang Lu's blog post: https://zhiganglu.com/post/pheatmap_change_annotation_colors/
    dimH <- dim(HeatDelta.Data)[1]
    dimW <- dim(HeatDelta.Data)[2]
    X <- pheatmap(HeatDelta.Data*HeatDelta.Value_Multiplier, col = rampCols, na_col="#CCCCCC",
                  scale="none", breaks=ignore.breaks*HeatDelta.Value_Multiplier,
                  legend_breaks = HeatDelta.Legend_Breaks*HeatDelta.Value_Multiplier, # Genomic signature breaks:  see legend of Table on Pg 191 of Karlin et al. 1998 or in text of Karlin and Mrazek 1997, the latter which states "We distinguish extremes of dinucleotide relative abundances as follows: extremely high, symbolically 111, r*XY $ 1.50; very high, 11, 1.30 # r*XY , 1.50; significantly high, 1, 1.23 # r*XY , 1.30; marginally high, (1), 1.20 # r*XY , 1.23; extremely low, 222, r*XY # 0.50; very low, 22, 0.50 , r*XY # 0.70; significantly low, 2, 0.70 , r*XY # 0.78; marginally low, (2), 0.78 , r*XY # 0.81."
                  legend_labels = HeatDelta.Legend_Labels,
                  cluster_rows = HeatDelta.Cluster_Data, annotation_row = HeatDelta.Group_Name, annotation_names_row=TRUE, annotation_colors = HeatDelta.Group_Colour,
                  cluster_cols = HeatDelta.Cluster_Data, annotation_col = HeatDelta.Group_Name, annotation_names_col=TRUE,
                  #cellwidth=20, cellheight=20,
                  display_numbers = TRUE, number_color = "black", fontsize_number = round(dimW/2.5), number_format ="%.0f",
                  fontsize_col = 10, main = paste0("k=",kmer, " absolute difference relative abundance * ", HeatDelta.Value_Multiplier), 
                  angle_col = "90"
    )
    print(X)
    HeatDelta.Width = tmp.W; HeatDelta.Height = tmp.H
    if(HeatDelta.Image_Output!=""){ dev.off() }
    
    rm(ignore.breaks, ignore.breakpoint, tmp.H, tmp.W)
    
  }
  
  stopTime = Sys.time()
  cat(paste0("\n\n___  HEATMAP (DELTA) RUNTIME INFO  ___\n\tRun start: ", startTime, " ", format(startTime,"%Z (UTC%z)"),
             "\n\tRun stop: ", stopTime, " ", format(stopTime,"%Z (UTC%z)"),
             "\n\tRun time: ", format(stopTime-startTime),"\n"), file = resultOutput.info)
  rm(startTime, stopTime)
  
} # Delta heatmap section


# Future to-do? Implement Chi^2 and G-test statistics check? (see: https://www.rdocumentation.org/packages/Biostrings/versions/2.40.2/topics/dinucleotideFrequencyTest)

#### FINAL SECTION ####
runtimeStop=Sys.time()
cat(paste0("\n\n___  OVERALL RUNTIME INFO  ___",
           "\n\tRun start: ", runtimeStart, " ", format(runtimeStart,"%Z (UTC%z)"),
           "\n\tRun stopped: ", runtimeStop, " ", format(runtimeStop,"%Z (UTC%z)"),
           "\n\tRun time: ",format(runtimeStop - runtimeStart),
           
           "\n\n\n___ SESSION INFORMATION ___\n\n"), file = resultOutput.info)

# Output sessionInfo() using sink() as learned through: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
sink(file = resultOutput.info, append = TRUE)
print(sessionInfo())
sink()
close(resultOutput.info)

# Stop running everything in the background
stopCluster(pClust); rm(pClust)  #  Stop the parallel cluster since we no longer need them now
closeAllConnections()  #  Safety catch in case any files did not close properly
cat("::: Finished running oNTORRA :::\n")

#stop("Finished!")
############  END OF MAIN SCRIPT ###################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####   DUMMY SET FREQ. DATA    #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# warning("WARNING: Dummy dataset [AATCCGGAT] is being used!")					
# # https://www.bioinformatics.nl/cgi-bin/emboss/compseq					
# ###  f(k=1) [ NOT f*!]					
#   #
#   # Output from 'compseq'
#   #
#   # The Expected frequencies are calculated from the observed single
#   # base or residue frequencies in these sequences
#   #
#   # The input sequences are:
#   #
# 
# 
# Word size	1
# Total count	9
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# A	3		0.3333333	0.3333333	1
# C	2		0.2222222	0.2222222	1
# G	2		0.2222222	0.2222222	1
# T	2		0.2222222	0.2222222	1
# 
# Other	0		0	0	10000000000
# 
# 
# #https://www.bioinformatics.nl/cgi-bin/emboss/compseq
# ###  f*(k=1)					
# #   #					
# #					
# # Output from 'compseq'					
# #					
# # The Expected frequencies are calculated from the observed single					
# # base or residue frequencies in these sequences					
# #					
# # The input sequences are:					
# #					
# 
# 
# Word size	1				
# Total count	18				
# 
# #					
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency	
# #					
# A	5		0.2777778	0.2777778	1
# C	4		0.2222222	0.2222222	1
# G	4		0.2222222	0.2222222	1
# T	5		0.2777778	0.2777778	1
# 
# Other	0		0	0	10000000000
# #					
# ###  f*(k=2)					
# #
# # Output from 'compseq'
# #
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	2
# Total count	16
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# AA	1	0.0625	0.0771605	0.81
# AC	0	0	0.0617284	0
# AG	0	0	0.0617284	0
# AT	4	0.25	0.0771605	3.24
# CA	0	0	0.0617284	0
# CC	2	0.125	0.0493827	2.53125
# CG	2	0.125	0.0493827	2.53125
# CT	0	0	0.0617284	0
# GA	2	0.125	0.0617284	2.025
# GC	0	0	0.0493827	0
# GG	2	0.125	0.0493827	2.53125
# GT	0	0	0.0617284	0
# TA	0	0	0.0771605	0
# TC	2	0.125	0.0617284	2.025
# TG	0	0	0.0617284	0
# TT	1	0.0625	0.0771605	0.81
# 
# Other	0		0	0	10000000000
# #					
# ###  f*(k=3), zero skipped					
# #
# # Output from 'compseq'
# #
# # Words with a frequency of zero are not reported.
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	3
# Total count	14
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# AAT	1		0.0714286	0.0214335	3.3325714
# ATC	2		0.1428571	0.0171468	8.3314286
# ATT	1		0.0714286	0.0214335	3.3325714
# CCG	2		0.1428571	0.0109739	13.0178571
# CGG	2		0.1428571	0.0109739	13.0178571
# GAT	2		0.1428571	0.0171468	8.3314286
# GGA	2		0.1428571	0.0137174	10.4142857
# TCC	2		0.1428571	0.0137174	10.4142857
# 
# Other	0		0	0	10000000000
# 
# 
# ###  f*(k=4), zero skipped					
# #
# # Output from 'compseq'
# #
# # Words with a frequency of zero are not reported.
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	4
# Total count	12
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# AATC	1		0.0833333	0.004763	17.496
# ATCC	2		0.1666667	0.0038104	43.74
# CCGG	2		0.1666667	0.0024387	68.34375
# CGGA	2		0.1666667	0.0030483	54.675
# GATT	1		0.0833333	0.004763	17.496
# GGAT	2		0.1666667	0.0038104	43.74
# TCCG	2		0.1666667	0.0030483	54.675
# 
# Other	0		0	0	10000000000
# 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Dummy set 2 :: CGATCCGGATA
# 
# 
# #
# # Output from 'compseq'
# #
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	1
# Total count	11
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# A	3		0.2727273	0.2727273	1
# C	3		0.2727273	0.2727273	1
# G	3		0.2727273	0.2727273	1
# T	2		0.1818182	0.1818182	1
# 
# Other	0		0	0	10000000000
# 
# # https://www.bioinformatics.nl/cgi-bin/emboss/compseq
# ##  f*(k=1), show all
# #
# # Output from 'compseq'
# #
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	1
# Total count	22
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# A	5		0.2272727	0.2272727	1
# C	6		0.2727273	0.2727273	1
# G	6		0.2727273	0.2727273	1
# T	5		0.2272727	0.2272727	1
# 
# Other	0		0	0	10000000000
# 
# ##  f*(k=2), zero skipped
# #
# # Output from 'compseq'
# #
# # Words with a frequency of zero are not reported.
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	2
# Total count	20
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# AT	4		0.2	0.0516529	3.872
# CC	2		0.1	0.0743802	1.3444444
# CG	4		0.2	0.0743802	2.6888889
# GA	3		0.15	0.0619835	2.42
# GG	2		0.1	0.0743802	1.3444444
# TA	2		0.1	0.0516529	1.936
# TC	3		0.15	0.0619835	2.42
# 
# Other	0		0	0	10000000000
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##  f*(k=3), zero skipped
# #
# # Output from 'compseq'
# #
# # Words with a frequency of zero are not reported.
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	3
# Total count	18
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# ATA	1		0.0555556	0.0117393	4.7324444
# ATC	3		0.1666667	0.0140872	11.8311111
# CCG	2		0.1111111	0.0202855	5.4773663
# CGA	1		0.0555556	0.0169046	3.2864198
# CGG	2		0.1111111	0.0202855	5.4773663
# GAT	3		0.1666667	0.0140872	11.8311111
# GGA	2		0.1111111	0.0169046	6.5728395
# TAT	1		0.0555556	0.0117393	4.7324444
# TCC	2		0.1111111	0.0169046	6.5728395
# TCG	1		0.0555556	0.0169046	3.2864198
# 
# Other	0		0	0	10000000000
# ##  f*(k=4), zero skipped
# 
# 
# Output file   outfile
# 
# #
# # Output from 'compseq'
# #
# # Words with a frequency of zero are not reported.
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	4
# Total count	16
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# ATCC	2		0.125	0.003842	32.5355556
# ATCG	1		0.0625	0.003842	16.2677778
# CCGG	2		0.125	0.0055324	22.5941358
# CGAT	1		0.0625	0.003842	16.2677778
# CGGA	2		0.125	0.0046103	27.112963
# GATA	1		0.0625	0.0032016	19.5213333
# GATC	2		0.125	0.003842	32.5355556
# GGAT	2		0.125	0.003842	32.5355556
# TATC	1		0.0625	0.0032016	19.5213333
# TCCG	2		0.125	0.0046103	27.112963
# 
# Other	0		0	0	10000000000
# 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# dummy dataset 3 [TGACGNATTAACNNNGGATC]					
# ##  f(k=1)  [NOT f*(k=1)]
# #
# #
# # Output from 'compseq'
# #
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	1
# Total count	20
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# A	5		0.25	0.25	1
# C	3		0.15	0.15	1
# G	4		0.2	0.2	1
# T	4		0.2	0.2	1
# 
# Other	4		0.2	0.2	1
# 
# 
# ##  f*(k=1)
# #
# # Output from 'compseq'
# #
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	1
# Total count	40
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# A	9		0.225	0.225	1
# C	7		0.175	0.175	1
# G	7		0.175	0.175	1
# T	9		0.225	0.225	1
# 
# Other	8		0.2	0.2	1
# 
# 
# ##  f*(k=2)
# #
# # Output from 'compseq'
# #
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	2
# Total count	38
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# AA	2		0.0526316	0.050625	1.0396361
# AC	2		0.0526316	0.039375	1.336675
# AG	0		0	0.039375	0
# AT	4		0.1052632	0.050625	2.0792723
# CA	1		0.0263158	0.039375	0.6683375
# CC	1		0.0263158	0.030625	0.8592911
# CG	2		0.0526316	0.030625	1.7185822
# CT	0		0	0.039375	0
# GA	3		0.0789474	0.039375	2.0050125
# GC	0		0	0.030625	0
# GG	1		0.0263158	0.030625	0.8592911
# GT	2		0.0526316	0.039375	1.336675
# TA	2		0.0526316	0.050625	1.0396361
# TC	3		0.0789474	0.039375	2.0050125
# TG	1		0.0263158	0.039375	0.6683375
# TT	2		0.0526316	0.050625	1.0396361
# 
# Other	12		0.3157895	0.2	1.5789474
# 
# 
# ##  f*(k=3), zeroes excluded
# #
# # Output from 'compseq'
# #
# # Words with a frequency of zero are not reported.
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	3
# Total count	36
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# AAC	1		0.0277778	0.0088594	3.1354105
# AAT	1		0.0277778	0.0113906	2.4386526
# ACG	1		0.0277778	0.0068906	4.0312421
# ATC	2		0.0555556	0.0088594	6.2708211
# ATT	1		0.0277778	0.0113906	2.4386526
# CGT	1		0.0277778	0.0068906	4.0312421
# GAC	1		0.0277778	0.0068906	4.0312421
# GAT	2		0.0555556	0.0088594	6.2708211
# GGA	1		0.0277778	0.0068906	4.0312421
# GTC	1		0.0277778	0.0068906	4.0312421
# GTT	1		0.0277778	0.0088594	3.1354105
# TAA	2		0.0555556	0.0113906	4.8773053
# TCA	1		0.0277778	0.0088594	3.1354105
# TCC	1		0.0277778	0.0068906	4.0312421
# TGA	1		0.0277778	0.0088594	3.1354105
# TTA	2		0.0555556	0.0113906	4.8773053
# 
# Other	16		0.4444444	0.2	2.2222222
# 
# 
# ##  f*(k=4), zeroes excluded
# #
# # Output from 'compseq'
# #
# # Words with a frequency of zero are not reported.
# # The Expected frequencies are calculated from the observed single
# # base or residue frequencies in these sequences
# #
# # The input sequences are:
# #
# 
# 
# Word size	4
# Total count	34
# 
# #
# # Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
# #
# ATCC	1		0.0294118	0.0015504	18.9705512
# ATTA	1		0.0294118	0.0025629	11.4760124
# CGTC	1		0.0294118	0.0012059	24.3907087
# GACG	1		0.0294118	0.0012059	24.3907087
# GATC	2		0.0588235	0.0015504	37.9411024
# GGAT	1		0.0294118	0.0015504	18.9705512
# GTCA	1		0.0294118	0.0015504	18.9705512
# GTTA	1		0.0294118	0.0019934	14.7548731
# TAAC	1		0.0294118	0.0019934	14.7548731
# TAAT	1		0.0294118	0.0025629	11.4760124
# TGAC	1		0.0294118	0.0015504	18.9705512
# TTAA	2		0.0588235	0.0025629	22.9520249
# 
# Other	20		0.5882353	0.2	2.9411765
