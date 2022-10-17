#!/bin/bash
# Obtain the nucleotide counts using compseq
# Moez Valliani 2020 NOV 03
#
# Note: This script does not have quotes around variable, ensure that your directory does NOT contain any spaces as it will break this or other scripts.
# (You can add quotes around variables, like "$inFile", to help "escape" spaces but it is best to avoid having spaces in general and use underscore instead).
# 
# USAGE:  01_*.bash  FASTA_FILE  [OUTPUT_DIRECTORY = oligonucleotide_count]
#
# This script essentially makes a temporary file containing the FASTA sequence (e.g. genome sequences), removes the header, then uses EMBOSS' compseq to obtain the nucleotide motif counts.
#
# Must have EMBOSS suite install and have compseq (apt install emboss). Used version 6.6.0?
# This script obtains the k=1 nucleotide frequency ONLY in the forward/sense direction. (To have both, define the initial doRev with "-reverse" instead of an empty string)
# The other kmers (k>1) obtains the nucleotide motifs in BOTH the sense and antisense directions


if [ "$(which compseq)" == "" ]; then echo "ERROR: Need to install EMBOSS compseq!">&2; exit; fi

# Safety check / input section
if [ $# -lt 1 ]; then echo "OPTIONS:  FASTA_FILE  [OUTPUT_DIRECTORY = oligonucleotide_count]">&2; exit; fi

inFile="$1"  # Reads the FASTA file used in the first argument/parameter
if [[ ! -f $inFile ]]; then echo "ERROR: Input file $inFile cannot be found!">&2; exit; fi

outDir="oligonucleotide_count"  # If there is a second argument/parameter, it will output to that directory, otherwise it will default to 'oligonucleotide_count'
if [ $# -ge 2 ]; then
	shift
	outDir="$*"
fi
if [[ ! -d $outDir ]]; then mkdir $outDir; fi

### Core


tmpFile="/tmp/oligonucleotide_count_$(date '+%y%m%d_%H%M%S')_$(basename $inFile)_$(md5sum $inFile | cut -f1 -d" ")"

# Pre-emptive modifications to the temp. file
# SED replaces all the FASTA headers in the file with a single >
# TR -d '\012'  removes all the newlines in the file
# TR '>' '\012' replaces all the single > with a new line (to separate the starting Kmer)
# sed -e "/^>/d" $inFile | tr -d '\012' > $tmpFile
# sed -e "s/^>.*/>/" $inFile | tr -d '\012' | tr '>' '\012' > $tmpFile

sed -e "s/^>.*//" $inFile > $tmpFile # | tr -d '\012' | sed -e "s/\s//g" > $tmpFile

doRev=""  # Used to keep track of when to reverse complement (when K > 1). If you want to have k=1 to account for both directions, set doRev="-reverse" instead of doRev=""
for k in {1..4}; do
	outFile=$outDir/k${k}${doRev}.$(basename $inFile).tsv
	compseq -sequence $tmpFile -word $k $doRev -calcfreq -outfile $outFile
	doRev="-reverse"  # Enables reverse complement count for k=2+
	
	# Check if other nucleotides are present in the NT frequency calculation, if so, append the counts to the end using awk
	if [ $k -eq 1 ] && [ $(tail -n1 $outFile | cut -f2) -ne 0 ]; then
		# Removes newlines and ATCG characters
		cat $tmpFile | tr -d '\012' | sed -e "s/[\sATCGatcg]//g" > $tmpFile.2
		awk '
			{
				$0=toupper($0)
				while($0){
					X = substr($0, 1, 1)    # Get the first character
					print X "\t" gsub(X,"") # gsub() returns the number of replaced characters
				}
			}
		' $tmpFile.2 >> $outFile
	fi
done

rm ${tmpFile}*