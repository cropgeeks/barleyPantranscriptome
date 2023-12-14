#!/bin/bash

#SBATCH --cpus-per-task=12
#SBATCH -o slurm-%x_%A.out 
#SBATCH --mem=80G

outputDir=`pwd`
sjFile=allSJs.filtered.tab
genomeFastaFile=Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.split.fasta

echo "generate index using SJ file"

STAR \
--runThreadN $SLURM_CPUS_PER_TASK \
--runMode genomeGenerate \
--genomeDir $outputDir \
--genomeFastaFiles $genomeFastaFile \
--sjdbFileChrStartEnd $sjFile

echo "done"

