#!/bin/bash


#script for variant calling with Freebayes

input_bam=$1
base_name=$2
ref_fasta=$3

#######################

workingDir=`pwd`
echo "workingDir = $workingDir"
echo "input_bam = $input_bam"
echo "TMPDIR = $TMPDIR"

echo "copy the input data to the local node"
cp $input_bam $TMPDIR
baiFile=`echo $input_bam | sed 's/bam/bai/g'`
cp $baiFile $TMPDIR
localBamFile=`basename $input_bam`
echo "localBamFile = $localBamFile"
cd $TMPDIR

#run the variant calling
source activate freebayes
echo "start Freebayes run"
date

freebayes \
-f $ref_fasta \
--min-alternate-count 3 \
--min-alternate-fraction 0.3 \
--ploidy 2 \
-m 0 \
-v $base_name.vcf \
--legacy-gls \
--use-duplicate-reads \
$localBamFile

echo "Freebayes run complete"
echo "local files in TMPDIR:"
ls -l

echo "copy the output data back to shared storage"
cp $base_name.vcf $workingDir

echo "done"
date


