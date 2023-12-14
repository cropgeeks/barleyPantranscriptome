#!/bin/bash

#SBATCH --cpus-per-task=32
#SBATCH -o slurm-%x_%j.out
#SBATCH --mem=120G
#SBATCH --partition=long

fileList=localFiles.txt
outputFile=barleyPanTranscriptome_variants_freebayes_splitcoords.bcf

date
echo "merge the VCF files"
bcftools merge \
--missing-to-ref \
--force-samples \
--merge all \
--threads $SLURM_CPUS_PER_TASK \
--file-list $fileList \
--output-type b \
--output $outputFile

echo "index BCF file"
bcftools index $outputFile

echo "done"
date
