#!/bin/bash

#SBATCH --mem=16G
#SBATCH -o slurm-%x_%A.out 


wholeChromFASTA=Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.fasta

splitChromFASTA=Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.split.fasta

echo "splitting pseudomolecules"

java \
-Xmx15g \
utils.fasta.SplitPseudomolecules \
$wholeChromFASTA \
$splitChromFASTA

echo "done"