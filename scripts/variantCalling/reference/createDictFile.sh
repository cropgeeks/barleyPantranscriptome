#!/bin/bash

#SBATCH --mem=5G

refseq=Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.split.fasta

java \
-jar $APPS/picard/picard.jar \
CreateSequenceDictionary \
R=$refseq \
O=$refseq.dict