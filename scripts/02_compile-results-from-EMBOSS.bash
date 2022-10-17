#!/bin/bash
# Compile the frequency results in to a single file
# Moez Valliani 2020 OCT 14
#
# USAGE: bash 02_*.bash
# 
#
# Note: This script only accounts for k=1-4. If you want to adjust, search for for k in {1..4} and change the range value

filePrefix="bac3"									  # If any prefix is appended to the filename
resultDir="$PWD/bac/1_nt-counts/*/*/oligonucleotide_count"   # The result directory containing the compseq results
outDir="$PWD/bac/2_oligocounts"	  # Output directory

if [[ $(ls $resultDir/k*.tsv 2>/dev/null | wc -l) -eq 0 ]]; then echo "ERROR: Cannot find the result directory containing the parsed k-mer files!">&2; exit; fi
if [[ ! -d "$outDir" ]]; then mkdir $outDir; fi

## NT frequency data
echo "Obtaining NT freq data"

echo "Obtaining oligomer freq data"

for k in {1..4}; do
	awk -v K=$k -v filePrefix=$filePrefix -F"\t" '
		# Oligo = The Oligoay containing the permutated values
		function iterATCG(Len, Pos, Str){
			# Len = oligo length
			# Pos = current position in the string
			# Str = the current string
			if(Pos > Len){
				Oligo[Str] = Str
				return
			}
		
			ACGT[1]="A";ACGT[2]="C";ACGT[3]="G";ACGT[4]="T";
			for(X in ACGT){
				iterATCG(Len, Pos+1, Str ACGT[X])
			}
		}
	
		FNR==1{split(FILENAME, tmp, "/"); ORG=tmp[length(tmp)-2]; ON[ORG]=ORG; next}  # Get the organism name when a new file is loaded based on the file path
		FNR<17{next}	# Skip the first 16 lines which are header info from compseq
		NF<2{next}		# Skip entries that have a blank line, which would occur after the entries have been reported
		{
			OG[ORG][$1]=$2     # OG[organism][k_motif]=count
			if($1=="Other"){nextfile}
		}

		END{
			# Generate the various ATCG oligomer combinations
			Header="Organism"
			
			# Obtaining all possible oligomer combination (4^k) for k = 2 to 4 (sets combinations to variable Oligo)
			iterATCG(K,1,"")


			# Making the oligomer count row result for each organism
			asorti(Oligo)
			for(i in Oligo){
				O=Oligo[i]
				Header=Header "\t" O

				for(ORG in OG){
					i=0
					if(O in OG[ORG]){ i=OG[ORG][O] }
					Res[ORG]=Res[ORG] "\t" i
				}
			}

			Header=Header "\tOther"
			for(ORG in OG){
				Res[ORG]=Res[ORG] "\t" OG[ORG]["Other"]
			}



			# Outputting the final combined results
			print Header

			asorti(ON)
			for(O in ON){
				ORG=ON[O]
				print ORG Res[ORG]
			}
		}
	' $resultDir/k${k}*.tsv > "$outDir/${filePrefix}-k${k}-oligo-freq.tsv"

done
