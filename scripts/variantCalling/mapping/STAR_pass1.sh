#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH -o slurm-%x_%A_%a.out 
#SBATCH --mem=60G
#SBATCH --partition=long
#SBATCH --array=1-22%11

#script for mapping to a reference using STAR
#@Micha Bayer, James Hutton Institute, 12 Feb 2021

#Since we are processing multiple samples we want to do the 2-pass mapping in accordance with the STAR recommendations
#This means doing the first pass mapping on all samples and then combining the splice junctions from all samples so they can be used as annotation for the second pass
#This will produce the best sensitivity for novel junctions.
#See these threads about this:
#https://groups.google.com/g/rna-star/c/VTX9TfapSfQ
#https://groups.google.com/g/rna-star/c/9C3W_BMfGXM/m/-rg7C6HURHsJ
#This is the script for pass1 - pass2 is in script STAR_pass2.sh


###############################VARIABLES###############################

#an array with all the sample names
samples[1]=Akashinriki
samples[2]=Barke
samples[3]=Chi_Ba_Damai
samples[4]=Clipper
samples[5]=Du_Li_Huang
samples[6]=FT11
samples[7]=Golden_Promise
samples[8]=Hockett
samples[9]=HOR_10350
samples[10]=HOR_13821
samples[11]=HOR_13942
samples[12]=HOR_21599
samples[13]=HOR_3081
samples[14]=HOR_3365
samples[15]=HOR_7552
samples[16]=HOR_8148
samples[17]=HOR_9043
samples[18]=Igri
samples[19]=Morex
samples[20]=OUN333
samples[21]=RGT_Planet
samples[22]=Stirling

#the raw data lives here
rawDataDir=rawData
#the name of the sample for this instance of the script
sampleName=${samples[SLURM_ARRAY_TASK_ID]}
#the path to this sample's data
sampleDir=$rawDataDir/$sampleName
#the dir with the index files for STAR - this uses Morex v3 (2021 long read based assembly)
starIndex=refseq/STAR_index


###############################GO PROCESS###############################
echo -e "\n================================="
echo "start processing sample $sampleName"
echo `date`

mkdir $sampleName
cd $sampleName

#STAR has to have a large numbers of files open during BAM sorting
#so we need to increase the ulimit -n value 
ulimit -n 10000

#these are the different tissue types we are dealing with
tissues=( Ca Co In Ro Sh )
#iterate over these
for tissue in ${tissues[*]}
do	
	echo -e "\n\n--------------process tissue $tissue"
	
	#iterate over the three bioreps for this tissue 
	#we want these all mapped separately so we can use them as internal QC control
	for biorep in {1..3}
	do
		echo -e "\n--------------process biorep $biorep for tissue $tissue"
		
		#prefix for our output
		outFileNamePrefix=$sampleName"_"$tissue""$biorep"_pass1_"
		
		#local temp dir for faster writing of output files; this gets deleted after the job ends, so we need to move the results to somewhere permanent before the job exits
		output_dir=$TMPDIR/$outFileNamePrefix
		echo "output_dir = $output_dir"
		mkdir $output_dir
		ls -lh $output_dir
		
		#the final destination for the output that we copy from the temp dir
		final_output_dir=`pwd`
		echo "final_output_dir = $final_output_dir"		
			
		#raw file names look like this:
		#Hor10350_Sh2_1.fq.gz
		#Hor10350_Sh2_2.fq.gz
		infileR1=`ls -1 $sampleDir/*$tissue""$biorep"_1.fq.gz"`
		infileR2=`ls -1 $sampleDir/*$tissue""$biorep"_2.fq.gz"`
		echo "input files:"
		echo $infileR1
		echo $infileR2

		echo "run STAR mapper"
		STAR \
		--genomeDir $starIndex \
		--readFilesIn $infileR1 $infileR2 \
		--runThreadN $SLURM_CPUS_PER_TASK \
		--outBAMsortingThreadN $SLURM_CPUS_PER_TASK \
		--outFileNamePrefix $output_dir/$outFileNamePrefix \
		--outBAMcompression 10 \
		--outSAMattrRGline ID:$sampleName""$tissue""$biorep SM:$sampleName LB:$tissue""$biorep \
		--twopassMode None \
		--alignIntronMin 60 \
		--alignIntronMax 15000 \
		--alignMatesGapMax 2000 \
		--alignEndsType Local \
		--alignSoftClipAtReferenceEnds No \
		--outSAMprimaryFlag AllBestScore \
		--outFilterMismatchNoverLmax 0.02 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 1 \
		--outFilterMatchNmin 0 \
		--outFilterMatchNminOverLread 0 \
		--outFilterMultimapNmax 15 \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate \
		--alignTranscriptsPerReadNmax 30000 \
		--readFilesCommand zcat \
		--outReadsUnmapped Fastx \
		--alignSJoverhangMin 7 \
		--alignSJDBoverhangMin 7 \
		--alignSJstitchMismatchNmax 0 1 0 0
		
		echo "moving output dir $output_dir to permanent storage in $final_output_dir"
		mv $output_dir $final_output_dir
	
	done

done

echo -e "\ndone processing sample $sampleName"
echo `date`

###############################END###############################

