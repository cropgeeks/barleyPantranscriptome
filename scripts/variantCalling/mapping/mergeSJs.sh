#!/bin/bash

#SBATCH -o slurm-%x_%A.out 
#SBATCH --mem=500M

#script for merging the splice junctions from a set of samples following STAR mapping pass 1
#this removes redundancy of splice junctions and filters for read support of at least 2 reads
#@Micha Bayer, James Hutton Institute, 19 Feb 2021

#from the STAR manual:
# SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. Note that STAR defines the junction start/end as intronic bases, while many other software define them as exonic bases. The columns have the following meaning:
# •	column 1: chromosome
# •	column 2: first base of the intron (1-based)
# •	column 3: last base of the intron (1-based)
# •	column 4: strand (0: undefined, 1: +, 2: -)
# •	column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
# •	column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
# •	column 7: number of uniquely mapping reads crossing the junction
# •	column 8: number of multi-mapping reads crossing the junction
# •	column 9: maximum spliced alignment overhang

echo "merge splice junctions from all samples"

cat */*/*pass1_SJ.out.tab \
| awk '($5 > 0 && $7 > 2 && $6==0)' \
| cut -f1-6 \
| sort \
| uniq \
> allSJs.filtered.tab

echo "done"