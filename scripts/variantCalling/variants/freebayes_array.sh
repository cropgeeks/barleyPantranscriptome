#!/bin/bash

#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH --mem=16G
#SBATCH --array=1-326%50
#SBATCH --partition=long

#a config file with all our BAM files listed
configFile=allBAMfiles.txt

#the reference sequence used for mapping
refSeq=Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.split.fasta


#parse an input config file with file paths to all our BAM files
#split the input by line
IFS=$'\n'
#read the whole file with the sample info into an array
lines=( `cat "$configFile" `)
#based on the SLURM_ARRAY_TASK_ID we pick a single BAM file from the input file, and this will be the one that gets used
#note the SLURM_ARRAY_TASK_ID starts at 1, but the lines array is zero-based, hence need to subtract 1 
bamFile=${lines[$SLURM_ARRAY_TASK_ID-1]}

#extract the name of the line and tissue for labeling of our output
sample=`basename $bamFile | sed 's/.recalibrated.bam//g' `

#run the script that implements the GATK small variant calling pipeline for RNAseq
bash \
./freebayes.sh \
$bamFile \
$sample \
$refSeq