---
title: "PanBaRT20"
author: "Wenbin Guo"
date: "2023-11-07"
output: html_document
---


# PanBaRT20 analysis pipeline

<hr>

Table of contents
-----------------

-   [Description](#description)
-   [Short-read RTD construction](#short-read-rtd-construction)
    +   [RNA-seq data pre-processing](#rna-seq-data-pre-processing)
    +   [RNA-seq read mapping](#rna-seq-read-mapping)
    +   [RNA-seq transcript assembly](#rna-seq-transcript-assembly)
    +   [Merge transcripts across samples](#merge-transcripts-across-samples)
-   [Long-read RTD construction](#long-read-rtd-construction)
    +   [Iso-seq data pre-processing](#iso-seq-data-pre-processing)
    +   [Iso-seq read mapping](#iso-seq-read-mapping)
    +   [Iso-seq transcript assembly](#iso-seq-transcript-assembly)
    +   [Assembly quality control](#assembly-quality-control)
-   [Genotype-specific RTDs](#genotype-specific-rtds)
-   [Barley PanBaRT20](#barley-panbart20)
    +   [Construction of linear pan-genome](#construction-of-linear-pan-genome)
    +   [Map GsRTDs to linear pan-genome](#map-gsrtds-to-linear-pan-genome)
    +   [Gene and trasncript association between PanBaRT20 and GsRTDs](#gene-and-trasncript-association-between-panbart20-and-gsrtds)
-   [References](#references)

<div align="justify">



Description
-----------


Short-read RTD construction
-----------

### RNA-seq data pre-processing

Fastp v0.20.1 (Chen et al., 2018) was used to remove adapters and filter reads with a quality score of less than 20 and a length of less than 30 bases.

```

# read1: read 1 of the paired-end RNA-seq reads
# read2: read 2 of the paired-end RNA-seq reads
# output_folder: the output directory to save the results

mkdir -p $output_folder

fastp \
-i ${read1} \
-I ${read2} \
-o ${output_folder}/trimmed_$(basename $read1) \
-O ${output_folder}/trimmed_$(basename $read2) \
-q 20 \
--cut_front \
--cut_tail \
-l 30 \
-h ${output_folder}/${sample}_fastp.html \
-j ${output_folder}/${sample}_fastp.json

```
<a href='#table-of-contents'>Go back to Table of contents</a>

### RNA-seq read mapping

The trimmed RNA-seq reads were mapped to the reference genome by using STAR v2.7.8a (Dobin et al., 2013; Dobin and Gingeras, 2015).

#### STAR mapping index

```

# genome_fasta: genome fasta file
# output_dir: output directory of the STAR index
# sj_file: The splice junction file for the pass 2 index. 

###---> index star mapping pass1
STAR \
--runMode genomeGenerate \
--genomeDir $output_dir \
--genomeFastaFiles $genome_fasta \
--outFileNamePrefix $output_dir \
--limitGenomeGenerateRAM 240000000000


###---> index star mapping pass2
STAR \
--runMode genomeGenerate \
--genomeDir $output_dir \
--genomeFastaFiles $genome_fasta \
--outFileNamePrefix $output_dir \
--limitGenomeGenerateRAM 240000000000 \
--sjdbFileChrStartEnd $sj_file

```

#### STAR mapping

```

###---> mapping scripts for pass1 and pass2

# Index_dir: STAR index directory
# read1, read2: read 1 and read 2 of the trimmed RNA-seq reads

STAR \
--genomeDir $Index_dir \
--readFilesIn $read1 $read1 \
--sjdbOverhang 100 \
--alignIntronMin 60 \
--alignIntronMax 15000 \
--alignMatesGapMax 2000 \
--alignEndsType Local \
--alignSoftClipAtReferenceEnds Yes \
--outSAMprimaryFlag AllBestScore \
--outFilterType BySJout \
--outFilterMismatchNmax 0 \
--outFilterMismatchNoverLmax 0.3 \
--outFilterScoreMinOverLread 0.66 \
--outFilterMatchNmin 0 \
--outFilterScoreMin 0 \
--outFilterMultimapNmax 15 \
--outFilterIntronMotifs RemoveNoncanonical \
--outFilterIntronStrands RemoveInconsistentStrands \
--outSJfilterReads All \
--outSJfilterCountUniqueMin -1 5 5 5 \
--outSJfilterCountTotalMin -1 5 5 5 \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--alignTranscriptsPerReadNmax 30000 \
--twopassMode None \
--readFilesCommand zcat \
--outReadsUnmapped Fastx \
--outFileNamePrefix $outfolder \
--outTmpDir $outTmpDir \
--alignSJoverhangMin 5 \
--alignSJDBoverhangMin 3 \
--outSJfilterOverhangMin -1 12 12 12 \
--outFilterMatchNminOverLread 0.66 \
--outFilterMismatchNoverReadLmax 1 \
--alignSJstitchMismatchNmax 0 0 0 0


###---> merge sj files from star mapping pass1
# sj_file: the splice junctions were saved in this file for the second pass of STAR index

cat star_result_pass1/*/SJ.out.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > $sj_file

```
<a href='#table-of-contents'>Go back to Table of contents</a>

### Transcript assembly

For complementary benefits, two assemblers, Stringtie v2.1.5 (Pertea et al., 2015) and Scallop v0.10.5 (Shao and Kingsford, 2017), were used to assemble the transcripts in the samples. 

#### Stringtie

```

# bam_file: the list of bam files of the read alignment
# sample: the sample ID
# save_dir: the output directory of the assembly

stringtie $bam_file \
-o ${save_dir}Stringtie_${sample}.gtf \
--rf \
-a 10 \
-c 2.5 \
-f 0 \
-g 50 \
-j 0.1 \
-M 1

```

#### Scallop

```

# bam_file: the list of bam files of the read alignment
# sample_id: the sample ID
# save_dir: the output directory of the assembly

scallop \
-i $bam_file \
-o ${save_dir}Scallop_${sample_id}.gtf \
--library_type first

```
<a href='#table-of-contents'>Go back to Table of contents</a>

### Merge transcripts across samples
We used RTDmaker (https://github.com/anonconda/RTDmaker) to merge the assembled transcripts across samples in each genotype. 


```

# rtd_dir: the directory of the Stringtie and Scallop assembly results
# sj_dir: the splice junction outputs in the STAR mapping pass2
# genome_fasta: the reference genome fasta file
# read_dir: the directory of the trimmed RNA-seq reads
# cultivar: genotype ID as the prefix of output file

python /home/wguo/scratch/pantrans/code/RTDmaker/RTDmaker.py ShortReads \
--assemblies $rtd_dir \
--SJ-data $sj_dir \
--genome $genome_fasta \
--fastq $read_dir \
--SJ-reads 10 2 \
--tpm 1 2 \
--fragment-len 0.7 \
--antisense-len 0.5 \
--add uniform \
--ram 8 \
--outpath $TMPDIR \
--outname ${cultivar} \
--prefix $cultivar \
--keep intermediary


```

<a href='#table-of-contents'>Go back to Table of contents</a>

Long-read RTD construction
-----------

### Iso-seq data pre-processing
IsoSeqv3 pipeline (https://github.com/PacificBiosciences/IsoSeq) was used to perform data pre-processing of Iso-seq long read data.

#### Circular Consensus Sequence calling
Circular consensus sequence (CCS) reads were generated using ccs v6.3.0. 

```

# read: Iso-seq subread bam file
# output: output directory

ccs --min-rq 0.9 -j $read $output

```

#### Primer removal and demultiplexing
Barcode sequence and primer sequences were identified and removed using lima v2.5.0. Read sequences were oriented from 5' to 3'. 

```

# input_bam: the output bam file from ccs read step
# barcode_fasta: barcode fasta file (see https://github.com/wyguo/PanBaRT20/blob/main/data/primer_IsoExpressKit.fasta)
# output_bam: output bam file

lima ${input_bam} ${barcode_fasta} ${output_bam} \
--isoseq --peek-guess --dump-clips

```

#### Isoseq3 refine
Reads were further refined by trimming ploy(A) tails and removing concatemers using isoseq refine v3.4.0, resulting in full-length non-concatemer (FLNC) reads.

```
# input_bam: the output bam file from lima step
# barcode_fasta: barcode fasta file (see https://github.com/wyguo/PanBaRT20/blob/main/data/primer_IsoExpressKit.fasta)
# output_bam: output bam file

isoseq3 refine ${input_bam} ${barcode_fasta} ${output_bam} --require-polya

###---> convert bam files to fasta files
samtools fasta ${output_bam} > ${output_fasta}

```
<a href='#table-of-contents'>Go back to Table of contents</a>

### Iso-seq read mapping

FLNC reads were mapped to reference genome by using Minimap2 v2.24 (Li, 2018). The maximum intron size was set to 15,000. 

```

# genome: reference genome fasta file
# input: FLNC read
# output: output directory of the mapping

minimap2 -ax splice:hq -uf -G $maxintron ${genome} ${input} -o ${output}

```
<a href='#table-of-contents'>Go back to Table of contents</a>

### Iso-seq transcript assembly

TAMA collapse generated transcriptome models from the read mapping results (Kuo et al., 2020). 

```

# input: the read mapping result
# genome: reference genome
# output: output directory

python /home/wguo/scratch/isoseq/app/tama/tama_collapse.py \
-s ${input} \
-f ${genome} \
-p ${output} \
-d merge_dup -x capped -m 0 -a 0 -z 0 -sj sj_priority -lde 30 -sjt 30 -rm low_mem

```

<a href='#table-of-contents'>Go back to Table of contents</a>

### Assembly quality control
#### Splice junction analysis


In-house R scripts are used to perform splice junction (SJ) analysis, where the focus is on detecting and handling mapping errors around SJs. The mapping error information is extracted from the TAMA outputs, allowing for the identification of transcripts with SJs that exhibit mapping errors within a range of +/- 10 base pairs around the SJs. In most cases, these problematic transcripts are removed from further analysis.

However, a special case arises when an SJ exhibits mapping errors, but it has still been positively identified as an SJ in the short-read STAR mapping results. In this unique scenario, such an SJ is considered a true SJ, and the corresponding transcript is spared from removal in the analysis. 

```
# input_dir: The data directory of TAMA output
# output_dir: The directory to save results
# genome_fasta: The name of genome fasta sequence file. If the file is not in the working directory, please provide the full path.
# sj_overhang An integer of SJ overhang. A SJ is filtered if it has mismatches at its left (-sj_overhang) or right (+sj_overhang) exonic regions.
# sjdatabase: A tab separated file of splice junction (SJ) information. In this analysis, the SJ information is from the STAR output of short-read mapping. The SJ files of multiple samples were mereged into a single file. If a canonical SJ in the isoseq assembly has a match in the database, it will be kept even though it has mismatches in the overhangs.
# output_dir: Results are saved in output_dir
 
##############################################################
library("rtracklayer")
library("GenomicRanges")
library("Biostrings")
library("tidyr")

start.time <- Sys.time()
message('|==========> Splice junction analysis: ',Sys.time(),' <==========|')
message('Step1: Check alignment error +/-',sj_overhang,' bases around splice junctions')

if(is.null(output_dir))
  output_dir <- input_dir

if(!file.exists(output_dir))
  dir.create(path = output_dir,recursive = T)

##############################################################
###---> sj error analysis
file_error <- list.files(path = input_dir,
                         pattern = '_local_density_error.txt$',
                         full.names = T,
                         recursive = T)
samples <- names(file_error) <- gsub('_flnc_collapsed_local_density_error.txt',
                                     '',
                                     basename(file_error))


file_polya <- list.files(path = input_dir,
                         pattern = '_polya.txt$',
                         full.names = T,
                         recursive = T)
names(file_polya) <- gsub('_flnc_collapsed_polya.txt',
                          '',
                          basename(file_polya))

file_reads <- list.files(path = input_dir,
                         pattern = '_trans_read.bed$',
                         full.names = T,
                         recursive = T)
names(file_reads) <- gsub('_flnc_collapsed_trans_read.bed',
                          '',
                          basename(file_reads))

file_trans <- list.files(path = input_dir,
                         pattern = '_collapsed.bed$',
                         full.names = T,
                         recursive = T)
names(file_trans) <- gsub('_flnc_collapsed.bed',
                          '',
                          basename(file_trans))

input_files <- data.frame(samples=samples,
                          sj_error=file_error[samples],
                          polya=file_polya[samples],
                          reads=file_reads[samples],
                          trans=file_trans[samples],
                          row.names = samples)

if(nrow(input_files)==1){
  input_files$trans_merged <- input_files$trans
} else {
  file_trans_merged <- list.files(path = input_dir,
                                  pattern = 'sample_merged.bed$',
                                  full.names = T,
                                  recursive = T)
  
  if(length(file_trans_merged)==0 & nrow(input_files)>0)
    stop('It seems that the dataset has multiple samples. Please run tama merge to merge
       transcripts from multipele samples.')
  
  input_files$trans_merged <- file_trans_merged
}

write.csv(input_files,file = file.path(output_dir,'input_files.csv'),row.names = F)
save(input_files,file = file.path(output_dir,'input_files.RData'))

##################################################################
###---> get sj label
# s <- 's1'
sj_trans <- lapply(samples, function(s){
  file2read <- input_files[s,'trans']
  reads_bed <- import(file2read)
  # reads_bed$gene_id <- gsub(';.*','',reads_bed$name)
  # reads_bed$transcript_id <- gsub('.*;','',reads_bed$name)
  
  ## get exons of the transcriptome
  exons <- unlist(blocks(reads_bed))
  
  ## add gene and transcript ids
  exons$gene_id <- gsub(';.*','',names(exons))
  exons$transcript_id <- gsub('.*;','',names(exons))
  names(exons) <- NULL
  introns <- exon2intron(exons,sorted = F)
  label <- data.frame(gene_id=paste0(s,'_',introns$gene_id),
                      transcript_id=paste0(s,'_',introns$transcript_id),
                      intron_id=paste0(s,'_',introns$transcript_id,'.',introns$num_introns),
                      sj=introns$label)
  label[!duplicated(label),]
})
names(sj_trans) <- samples

##################################################################
###---> sj database
## sjdatabase
if(!is.null(sjdatabase)){
  if(!file.exists(sjdatabase))
    stop("The sj database file does not exist in the data direcotry.")
  message('  ->Reading SJ database: ',basename(sjdatabase),'\n')
  sj <- read.table(sjdatabase, header=FALSE,sep = '\t',quote = "\"",dec='.',fill = T,comment.char = "")
  sj <- sj[,1:4]
  sj <- sj[!duplicated(sj),]
  colnames(sj) <- c('seqnames','start','end','strand')
  sj$strand[sj$strand==0] <- '*'
  sj$strand[sj$strand==1] <- '+'
  sj$strand[sj$strand==2] <- '-'
  intron_db <- paste0(sj$seqnames,':',sj$start,'-',sj$end,';',sj$strand)
  intron_db <- unique(intron_db)
} else {
  intron_db <- NULL
}

##################################################################
###---> sj error
p <- paste0(c(rep('_',sj_overhang),'>',rep('_',sj_overhang)),collapse = '')

sj2keep <- list()
sj2remove <- list()
for(s in samples){
  file2read <- input_files[s,'sj_error']
  density_error <- read.table(file2read,header = T,sep="\t",quote = "\"",dec=".")
  density_error <- density_error[,c("cluster_id","scaff_name","start_pos","end_pos","strand","sj_error_simple" )]
  
  ## Filter read with "na" sj_error_simple
  idx <- which(density_error$sj_error_simple=='na')
  density_error <- density_error[-idx,]
  
  ## separate sj_error_simple to multiple rows with deliminate ";"
  density_error <- separate_rows(density_error, sj_error_simple,sep = ';')
  density_error$num_introns <- sequence(rle(density_error$cluster_id)$length)
  # density_error$sample <- s
  
  
  ###########
  ## only read column 4
  file2read <- input_files[s,'reads']
  colClasses <- rep("NULL",12)
  colClasses[4] <- "character"
  
  trans <- read.table(file2read, header=FALSE,sep="\t",quote = "\"",dec=".",
                      colClasses = colClasses)
  trans <- as.vector(t(trans))
  
  ## get cluster and transcript ids
  name.idx <- gsub(';.*','',trans)
  names(name.idx) <- gsub('.*;','',trans)
  density_error$transcript_id <- paste0(s,'_',name.idx[density_error$cluster_id])
  density_error$intron_id <- paste0(density_error$transcript_id,'.',density_error$num_introns)
  
  ## add sj labels
  sj_db <- sj_trans[[s]]
  idx <- match(density_error$intron_id,sj_db$intron_id)
  density_error$sj <- sj_db$sj[idx]
  
  #########
  idx1 <- grep(p,density_error$sj_error_simple)
  idx2 <- which(density_error$sj %in% intron_db)
  idx2keep <- union(idx1,idx2)
  kep <- unique(density_error$sj[idx2keep])
  rem <- setdiff(density_error$sj,kep)
  sj2keep <- c(sj2keep,setNames(object = list(kep),nm = s))
  sj2remove <- c(sj2remove,setNames(object = list(rem),nm = s))
}
sj_error <- setdiff(Reduce(union,sj2remove),Reduce(union,sj2keep))

message('Step2: Check nocanonical splice junctions')
###############################################################
###---> check canonical sj
# if(!file.exists(file2read))
#   stop('File "*_collasped.bed" is missing in the data directory:\n',data_dir)
file2read <- input_files$trans_merged[1]
trans_bed2filter <- import(file2read)
trans_bed2filter$gene_id <- gsub(';.*','',trans_bed2filter$name)
trans_bed2filter$transcript_id <- gsub('.*;','',trans_bed2filter$name)
save(trans_bed2filter,file=file.path(output_dir,'trans_bed2filter.RData'))

## get exons of the transcriptome
exons <- unlist(blocks(trans_bed2filter))

## add gene and transcript ids
exons$gene_id <- gsub(';.*','',names(exons))
exons$transcript_id <- gsub('.*;','',names(exons))
names(exons) <- NULL
introns <- exon2intron(exons,sorted = F)

## get the coordiantes of first two bases of introns
sj_start <- resize(introns, width=2, fix="start", use.names=T,ignore.strand=T)
sj_start <- split(sj_start,seqnames(sj_start))

## get the coordinates of the last two bases of introns
sj_end <- resize(introns, width=2, fix="end", use.names=T,ignore.strand=T)
sj_end <- split(sj_end,seqnames(sj_end))

message('  ->Reading reference genome: ',basename(genome_fasta),'\n')
###---> read the genome sequence fasta file
dna <- readBStringSet(filepath = genome_fasta)
names(dna) <- gsub(' .*','',names(dna))

chr <- levels(seqnames(exons))
names(chr) <- chr
dna <- dna[chr]

## extract sj motifs
sj_motif <- lapply(chr, function(x){
  # message(x)
  gr_start <- sj_start[[x]]
  if(is.null(gr_start))
    return(NULL)
  
  ## extract motif of first two bases of introns
  m_start <- extractAt(x = dna[[x]],at = gr_start@ranges)
  
  # message(x)
  gr_end <- sj_end[[x]]
  if(is.null(gr_end))
    return(NULL)
  
  ## extract motif of last two bases of introns
  m_end <- extractAt(x = dna[[x]],at = gr_end@ranges)
  
  ## put them into a data frame, add information of gene id, transcript id, strand
  motif <- data.frame(gene_id=gr_start$gene_id,
                      transcript_id=gr_start$transcript_id,
                      strand=strand(gr_start),
                      num_introns=gr_start$num_introns,
                      sj=gr_start$label,
                      motif=paste0(as.vector(m_start),as.vector(m_end)))
  motif
})

sj_motif <- do.call(rbind,sj_motif)
sj_motif$motif <- toupper(sj_motif$motif)
## transcript id + intron numbering as row names
# rownames(sj_motif) <- paste0(sj_motif$transcript_id,'.',sj_motif$num_introns)

idx <- which(sj_motif$strand=='-')
motif <- sj_motif[idx,'motif']
motif_unique <- unique(motif)

## reverse complement
motif_crev <- DNArevComplement(motif_unique)
names(motif_crev) <- motif_unique
sj_motif[idx,'motif'] <- motif_crev[motif]

# canonical <- c("GTAG","GCAG","ATAC","CTAC","CTGC","GTAT")
canonical <- c("GTAG","GCAG","ATAC")
sj_noncanonical <- unique(sj_motif$sj[!(sj_motif$motif %in% canonical)])

sj_filter <- list(sj_error=sj_error,sj_noncanonical=sj_noncanonical)
save(sj_filter,file=file.path(output_dir,'sj_filter.RData'))

message('Done: ',Sys.time())
end.time <- Sys.time()
time.taken <- end.time - start.time
message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')


######Defined function
#' Get introns from exon coordinates
#' @param exons A GRanges object of exons
#' @param tx2gene A data frame of transcript-gene association
#' @param add_label Logical, whether to add label of seqnames:start-end;strand
#' @return A GRanges object of introns
#'
exon2intron <- function(exons,tx2gene=NULL,add_label=T,sorted=T){
  
  ## filter gene with multiple strands
  strand_idx <- elementNROWS(range(GenomicRanges::split(exons,exons$gene_id)))
  strand_idx <- names(strand_idx)[strand_idx>1]
  exons <- exons[!(exons$gene_id %in% strand_idx),]
  
  if(is.null(tx2gene)){
    tx2gene <- DataFrame(transcript_id=exons$transcript_id,gene_id=exons$gene_id)
    tx2gene <- tx2gene[!duplicated(tx2gene),]
    rownames(tx2gene) <- tx2gene$transcript_id
  }
  # message('Split exons')
  exons <- GenomicRanges::split(exons,exons$transcript_id)
  # message('Get introns')
  introns <- psetdiff(unlist(range(exons),use.names=FALSE),exons)
  introns <- unlist(introns)
  if(NROW(introns)==0)
    return('No introns')
  
  introns$transcript_id <- names(introns)
  introns$gene_id <- tx2gene[introns$transcript_id,'gene_id']
  introns$feature <- 'intron'
  names(introns) <- NULL
  introns$num_introns <- sequence(rle(introns$transcript_id)$length)
  
  if(sorted)
    introns <- sort(introns,by=~seqnames + start + end)
  if(add_label)
    introns$label <- paste0(seqnames(introns),':',start(introns),'-',end(introns),';',strand(introns))
  return(introns)
}


intron2SJchain <- function(intron){
  gr <- GenomicRanges::split(intron,intron$transcript_id)
  s <- strand(gr)
  s <- unlist(runValue(s))
  
  chr <- seqnames(gr)
  chr <- unlist(runValue(chr))
  
  sj_chain <- paste0(start(intron),'-',end(intron))
  sj_chain <- GenomicRanges::split(sj_chain,intron$transcript_id)
  sj_chain <- sapply(sj_chain, function(x) paste0(x,collapse = ';'))
  chr <- chr[names(sj_chain)]
  s <- s[names(sj_chain)]
  label <- paste0(chr,':',sj_chain,':',s)
  names(label) <- names(sj_chain)
  label
}

getMultiExonTrans <- function(exon){
  gr <- GenomicRanges::split(exon,exon$transcript_id)
  n <- elementNROWS(gr)
  trans <- names(n)[n>1]
  exon[exon$transcript_id %in% trans,]
}

getMonoExonTrans <- function(exon){
  gr <- GenomicRanges::split(exon,exon$transcript_id)
  n <- elementNROWS(gr)
  trans <- names(n)[n==1]
  exon[exon$transcript_id %in% trans,]
}

getMonoExonGenes <- function(exon){
  gr <- GenomicRanges::split(exon,exon$gene_id)
  n <- elementNROWS(gr)
  genes <- names(n)[n==1]
  exon[exon$gene_id %in% genes,]
}

#getGeneLoci->getGeneRange
getGeneRange <- function(exon){
  gr <- GenomicRanges::split(exon,exon$gene_id)
  gr <- unlist(range(gr))
  gr$gene_id <- names(gr)
  names(gr) <- NULL
  # gr <- sort(gr,by=~seqnames+start+end)
  gr
}

#getTransCallapse->getGeneLoci
getGeneCollapse <- function(exon){
  gr <- GenomicRanges::split(exon,exon$gene_id)
  gr <- reduce(gr)
  gr <- unlist(gr)
  gr$gene_id <- names(gr)
  names(gr) <- NULL
  # gr <- sort(gr,by=~seqnames+start+end)
  gr
}

getOverlapRange <- function(exon){
  gr <- GenomicRanges::split(exon,exon$transcript_id)
  gr <- range(gr)
  gr <- unlist(gr)
  gr <- reduce(gr)
  gr$gene_id <- paste0(seqnames(gr),':',start(gr),';',strand(gr))
  gr
}


#getTransLoci->getTransRange
getTransRange <- function(exon){
  gr <- GenomicRanges::split(exon,exon$transcript_id)
  gr <- unlist(range(gr))
  gr$transcript_id <- names(gr)
  names(gr) <- NULL
  # gr <- sort(gr,by=~seqnames+start+end)
  gr
}

cleanMcols <- function(gr){
  meta <- mcols(gr)
  type_id <- ifelse('type' %in% colnames(meta),'type','feature')
  idx <- c(type_id,'transcript_id','gene_id')
  meta <- meta[,idx]
  mcols(gr) <- meta
  gr
}

subGR <- function(gr,cood){
  cood <- gsub(',','',cood)
  s <- unlist(strsplit(cood,split = ':|-'))
  gr[seqnames(gr)==s[1] & start(gr) > s[2] & end(gr) < s[3]]
}

```

<a href='#table-of-contents'>Go back to Table of contents</a>


#### Transcript start and end analysis

The sequencing frequency of transcription start sites (TSSs) and transcription end sites (TESs) is notably higher for highly expressed genes, with the number of reads supporting genuine TSSs/TESs being significantly distinguishable from random occurrences. In contrast, for lowly expressed genes, the reads are essentially indistinguishable from background noise when subjected to significance testing. As a result, we devised a two-step model to establish high-confidence TSS and TES sites for inclusion in the long-read transcriptome dataset.

Firstly, we employed a Binomial test to investigate the deviations between the observed reads at specific genomic sites (TSS/TES) and the expected reads based on randomness. For a gene with *s* genomic sites and a total of n mapped reads, the null hypothesis assumed that the read distribution was random. Each site had a fair probability (*x = 1/s*) of being sequenced, resulting in an expected number of reads for each site of *n/s*. We calculated the probability of observing *k* or more reads under this null hypothesis and defined a one-tailed test p-value as follows:

<math xmlns="http://www.w3.org/1998/Math/MathML" display="block" class="tml-display" style="display:block math;">
  <mrow>
    <mi>p</mi>
    <mo>=</mo>
    <mrow>
      <munderover>
        <mo movablelimits="false">∑</mo>
        <mrow>
          <mi>i</mi>
          <mo>=</mo>
          <mi>k</mi>
        </mrow>
        <mi>n</mi>
      </munderover>
    </mrow>
    <mrow>
      <mo fence="true">(</mo>
      <mfrac linethickness="0px">
        <mi>n</mi>
        <mi>i</mi>
      </mfrac>
      <mo fence="true">)</mo>
    </mrow>
    <msup>
      <mi>x</mi>
      <mi>i</mi>
    </msup>
    <mo form="prefix" stretchy="false">(</mo>
    <mn>1</mn>
    <mo>−</mo>
    <mi>x</mi>
    <msup>
      <mo form="postfix" stretchy="false">)</mo>
      <mrow>
        <mi>n</mi>
        <mo>−</mo>
        <mi>i</mi>
      </mrow>
    </msup>
  </mrow>
</math>

A significant p-value indicated the rejection of the null hypothesis that the genomic site occurred randomly. We adjusted these p-values using the Benjamini & Hochberg method to control the false discovery rate (FDR) associated with multiple hypotheses (Benjamini and Hochberg, 1995). To qualify as a high-confidence TSS/TES, it must meet the criteria of *k ≥ n/s* and *FDR < 0.05*. Recognizing the stochastic nature of TSS/TES, we considered reads whose ends were within ±50 nucleotides upstream and downstream of a significant genomic site as high-confidence.

Secondly, for genes without significant TSS/TES, a stringent threshold was applied, retaining genomic sites supported by a minimum of 2 reads within a sliding window of ±5 nucleotides. 

```
# Transcript start site (TSS) and end site (TES) analysis
# data_dir: The data directory of TAMA output.
# cut_off: A numeric value of BH adjsuted p-value (fdr), p-value or probablity cut-off.
# cut_at: The \code{cut_off} is applied to "fdr" (adjustd p-value), "pval" (p-value) or "prob" (probability).
# TSS_region,TES_region: An integer of enriched TSS/TES wobbling region. TSSs/TESs locate within enriched
# site+/-site_region are all treated as high confident sites. Default: 50.
# min_read: The minimum reads in the window around TSS/TES to support high confident sites. Default: 2.
# bin: An integer of window size to aggregate the reads upstream and downstream a TSS/TES.
# The total reads in [site-bin,site+bin] will be calculated. Default: 5.
# Results are saved in data_dir
# 

# cut_off = 0.05
# cut_at = 'fdr',
# TSS_region = 50
# TES_region = 50
# min_read = 2
# bin = 5


start.time <- Sys.time()
message('|=======> Transcript start and end site analysis: ',Sys.time(),' <=======|')
###########################################################################
### Step 1: Filter transcripts with polyA
###########################################################################
##---> get polya transcript
load(file.path(data_dir,'input_files.RData'))

message('Step 1: Get reads with potential poly A truncation')
### get read id of polya
cluster_polya=lapply(input_files$samples, function(s){
  ## input polya
  file2read <- input_files[s,'polya']
  polya <- read.table(file2read,
                      header=T,sep="\t",quote = "\"",dec=".")
  unique(polya$cluster_id)
})
cluster_polya <- unique(unlist(cluster_polya))

reads_bed <- lapply(input_files$samples, function(s){
  ## input the reads
  message(s)
  file2read <- input_files[s,'reads']
  bed <- import(file2read)
  bed$transcript_id <- paste0(s,'_',gsub(';.*','',bed$name))
  bed$gene_id <- gsub('[.].*','',bed$transcript_id)
  bed$cluster_id <- gsub('.*;','',bed$name)
  bed$name <- NULL
  bed
})
reads_bed <- do.call(c,reads_bed)
trans_reads <- table(reads_bed$transcript_id)
save(trans_reads,file = file.path(data_dir,'trans_reads.RData'))
###########################################################################
### Step 2: Generate statistics of TSS and TES
###########################################################################
message('\nStep 2: Transcript start site (TSS) analysis')
## generate TSS statistics
TSS_summary <- SummaryTranSite(gr = reads_bed,type = 'TSS',bin = bin)

message('\nStep 3: Transcript end site (TES) analysis')
## generate TES statistics
# transcript clusters with polya were removed
reads_bed_nopolya <- reads_bed[!(reads_bed$cluster_id %in% cluster_polya),]
TES_summary <- SummaryTranSite(gr = reads_bed_nopolya,type = 'TES',bin = bin)

save(TSS_summary,file = file.path(data_dir,'TSS_summary.RData'))
save(TES_summary,file = file.path(data_dir,'TES_summary.RData'))

###########################################################################
### Step 3: Extract high confident TSS and TES
###########################################################################
message('\nStep 4: Extract high confident TSS and TES')
site_stats <- TSS_summary
type <- ifelse('TSS' %in% colnames(site_stats),'TSS','TES')
HC_TSS <- HCsite(site_stats = site_stats,cut_off = cut_off,cut_at = cut_at,
                 site_region = TSS_region,min_read = min_read,type = type)

## Get HC TES
site_stats <- TES_summary
type <- ifelse('TSS' %in% colnames(site_stats),'TSS','TES')
HC_TES <- HCsite(site_stats = site_stats,cut_off = cut_off,cut_at = cut_at,
                 site_region = TES_region,min_read = min_read,type = type)

save(HC_TSS,file = file.path(data_dir,'HC_TSS.RData'))
save(HC_TES,file = file.path(data_dir,'HC_TES.RData'))

##########################################################################
end.time <- Sys.time()
time.taken <- end.time - start.time
message('\nDone: ',Sys.time())
message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')


######Defined function
#' Extract high confident transcript start or end site
#' @param site_stats A DataFrame of TSS or TES statistics generated by \code{SummaryTranSite}.
#' @param cut_off A numeric value of BH adjsuted p-value (fdr), p-value or probablity cut-off.
#' @param cut_at The \code{cut_off} is applied to "fdr" (adjustd p-value), "pval" (p-value) or "prob" (probability).
#' @param site_region An integer of enriched site region. TSSs/TESs locate within enriched
#' site+/-site_region are all treated as high confident sites.
#' @param min_read The minimum reads in a bin to support high confident sites.
#' @param type "TSS" for transcript start site and "TES" for transcript end site.
#' @return A DataFrame of enriched TSS or TES sites.
#'
HCsite <- function(site_stats,
                   cut_off=0.05,
                   cut_at=c('fdr','pval','prob'),
                   site_region=50,
                   min_read=2,
                   type=c('TSS','TES')){
  type <- match.arg(arg = type,choices = c('TSS','TES'))
  cut_at <- match.arg(arg = cut_at,choices = c('fdr','pval','prob'))

  if(!(type %in% colnames(site_stats)))
    stop('Please provide valid type: "TSS" or "TES"')

  #############################################################
  ###---> high abundance
  site_stats_enriched <- site_stats[site_stats[,cut_at] < cut_off &
                                      site_stats$read > site_stats$gene_read_mean,]
  site_stats_enriched$label2 <- paste0(site_stats_enriched$seqnames,';',site_stats_enriched$strand)
  site_stats_enriched_split <- GenomicRanges::split(site_stats_enriched,site_stats_enriched$label2)

  label_enriched <- lapply(names(site_stats_enriched_split), function(i){
    # print(i)
    site <- site_stats_enriched_split[[i]][,type]
    n <- rep((2*site_region+1),length(site))
    insite_region <- sequence(n,from = pmax(rep(0,length(site)),site-site_region))
    # insite_region <- insite_region[insite_region %in% site_stats[,type]]
    label <- paste0(i,';',insite_region)
    label <- label[which(label %in% site_stats$label)]
    label
  })
  label_enriched <- unique(unlist(label_enriched))
  site_high <- site_stats[site_stats$label %in% label_enriched,]

  #############################################################
  ###---> low abundance

  genes <- setdiff(site_stats$gene_id,site_high$gene_id)
  site_low <- site_stats[site_stats$gene_id %in% genes,]
  site_low <- site_low[site_low$read_in_bin >= min_read,]
  
  site_high$enriched <-'Yes'
  site_low$enriched <- 'No'
  
  result <- rbind(site_high,site_low)
  result <- sort(result,by=as.formula(paste0('~seqnames+',type)))
  result
}

```
<a href='#table-of-contents'>Go back to Table of contents</a>


#### Filter transcripts with low-quality SJs or TSS/TES 
Based on the splice junction and transcript start (TSS) and end (TES) analysis in previous steps, the transcripts with problematic or low-quality SJs or TSS/TES were filtered. 

```
# Filter collasped transcript datasets
# data_dir Data directory
# Results are saved in data_dir


message('|==========> Apply SJ and TSS/TES filters: ',Sys.time(),' <==========|')
message('Step 1: Load intermediate datasets\n')
start.time <- Sys.time()

if(!file.exists(file.path(data_dir,'HC_TSS.RData')) |
   !file.exists(file.path(data_dir,'HC_TES.RData')) |
   !file.exists(file.path(data_dir,'trans_bed2filter.RData')) |
   !file.exists(file.path(data_dir,'sj_filter.RData')) |
   !file.exists(file.path(data_dir,'input_files.RData'))
)
  stop('SJ and TSS/TES analysis must be conducted before RTD filter')

load(file.path(data_dir,'HC_TSS.RData'))
load(file.path(data_dir,'HC_TES.RData'))
load(file.path(data_dir,'trans_bed2filter.RData'))
load(file.path(data_dir,'sj_filter.RData'))
load(file.path(data_dir,'input_files.RData'))

message('Step 2: Apply splice junction filter\n')
##############################################################
###---> filter transcripts with error splice junction
sj_filter0 <- sj_filter
sj_filter <- unique(unlist(sj_filter))

exons <- unlist(blocks(trans_bed2filter))
## add gene and transcript ids
exons$gene_id <- gsub(';.*','',names(exons))
exons$transcript_id <- gsub('.*;','',names(exons))
names(exons) <- NULL
introns <- exon2intron(exons,sorted = F)
trans2filter <- unique(introns$transcript_id[introns$label %in% sj_filter])
gr <- trans_bed2filter[!(trans_bed2filter$transcript_id %in% trans2filter),]

trans_sj_error <- unique(introns$transcript_id[introns$label %in% sj_filter0$sj_error])
gr_sj_error <- trans_bed2filter[trans_bed2filter$transcript_id %in% trans_sj_error,]
genes_sj_error <- unique(gsub(';.*','',gr_sj_error$name))
trans_sj_error <- unique(gsub('.*;','',gr_sj_error$name))

trans_sj_noncanonical <- unique(introns$transcript_id[introns$label %in% sj_filter0$sj_noncanonical])
gr_sj_noncanonical <- trans_bed2filter[trans_bed2filter$transcript_id %in% trans_sj_noncanonical,]
genes_sj_noncanonical <- unique(gsub(';.*','',gr_sj_noncanonical$name))
trans_sj_noncanonical <- unique(gsub('.*;','',gr_sj_noncanonical$name))

message('Step 3: Apply TSS and TES filter\n')
##############################################################
###---> filter transcripts with low quality TSS and TES
TES_label <- unique(HC_TES$label)
TSS_label <- unique(HC_TSS$label)


## Match labels of TSS
gr_start <- resize(gr,fix = 'start',width = 1,ignore.strand=F)
gr_start_label <- paste0(seqnames(gr_start),';',strand(gr_start),';',start(gr_start))
trans_start <- gr_start$transcript_id[gr_start_label %in% TSS_label]

## Match labels of TES
gr_end <- resize(gr,fix = 'end',width = 1,ignore.strand=F)
gr_end_label <- paste0(seqnames(gr_end),';',strand(gr_end),';',end(gr_end))
trans_end <- gr_end$transcript_id[gr_end_label %in% TES_label]

## Subset rtd with HC transcripts
trans <- intersect(trans_start,trans_end)
rtd_tss_tes_filtered <- gr[!(gr$transcript_id %in% trans),]
trans_tss_tes_filtered <- unique(rtd_tss_tes_filtered$transcript_id)
genes_tss_tes_filtered <- unique(rtd_tss_tes_filtered$gene_id)

rtd_bed <- gr[gr$transcript_id %in% trans,]

stat_filterd <- data.frame(Category=c('SJ error',
                                      'SJ noncanonical',
                                      'TSS/TES filter'),
                           Genes=c(length(genes_sj_error),
                                   length(genes_sj_noncanonical),
                                   length(genes_tss_tes_filtered)),
                           Transcripts=c(length(trans_sj_error),
                                         length(trans_sj_noncanonical),
                                         length(trans_tss_tes_filtered)))
write.csv(stat_filterd,file=file.path(data_dir,'stat_filterd.csv'),row.names = F)

message('Step 4: Save the filtered RTD\n')
file2save <- input_files$trans_merged[1]

prefix <- paste0(gsub('.bed','',basename(file2save)),'_filtered')
export(object = gr_sj_error,file.path(data_dir,'transcript_sj_error.bed'))
export(object = gr_sj_noncanonical,file.path(data_dir,'transcript_sj_noncanonical.bed'))
export(object = rtd_tss_tes_filtered,
       file.path(data_dir,'transcript_tss_tes_filtered.bed'))

export(object = rtd_bed,file.path(data_dir,paste0(prefix,'.bed')),format = 'bed')

file2merge <- data.frame(file_name=file.path(normalizePath(data_dir,winslash = '/'),
                                             paste0(prefix,'.bed')),
                         cap_flag='capped',
                         merge_priority='1,1,1',
                         source_name=prefix)

write.table(file2merge, file = file.path(data_dir,'tama_merge_final_rtd_file.txt'),
            sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

##########################################################################

end.time <- Sys.time()
time.taken <- end.time - start.time
message('Done: ',Sys.time())
message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')

```

#### TAMA merge of redundant transcripts

To mitigate redundant transcripts resulting from stochasticity, transcripts with small variances within 50 nucleotides at the 5' and 3' UTR regions were removed using the TAMA merge tool with the parameters "-m 0 -a 50 -z 50 -d merge_dup."

```
# data_idr: the working directory of the long read assembly
# data_idr


python /home/wguo/scratch/isoseq/app/tama/tama_merge.py \
-f $data_idr/tama_merge_final_rtd_file.txt \
-p $data_idr/htsm_rtd \
-m 0 \
-a $t5 \
-z $t3 \
-d merge_dup
```
<a href='#table-of-contents'>Go back to Table of contents</a>

Genotype-specific RTDs
----------- 

The short-read assembly and long-read assembly in each genotype were merged to generate the genotype-specific RTDs (GsRTDs).

- The gene models of long-read RTD was used as the reference.
- The transcripts in short-read RTD that contribute to novel SJs or novel gene loci were integrated into the long-read RTD.
- In the integration result, the overlapped transcripts were assigned with the same gene ID. 
- If a set of overlapped transcripts can be divided into two groups and the maximum overlap length of these two group is < 5%, they are treated as separate gene models. 

```
# inf_file: short-read RTD gtf or bed file
# ref_file: long-read RTD gtf or bed file
# prefix_ref='ref': prefix attached to long-read RTD gene ids
# prefix_inf='inf': prefix attached to short-read RTD gene ids
# prefix_gene: NULL
# chimeric_tolerance: percentage cut-off of chimeric gene models
# genome_file: reference genome fasta file
# check_canonical: logical, whether to check canoncial SJs
# data_dir: the working direcotry


start.time <- Sys.time()
if(is.null(data_dir))
  data_dir <- getwd()

if(!file.exists(data_dir))
  dir.create(data_dir,recursive = T)

message('Step 1: Load gtf files')
## load gtf files
ref <- import(ref_file)
inf <- import(inf_file)

if(grepl(pattern = '[.]bed',inf_file)){
  inf <- unlist(blocks(inf))
  inf$gene_id <- gsub(';.*','',names(inf))
  inf$transcript_id <- gsub('.*;','',names(inf))
  inf$type <- 'exon'
  names(inf) <- NULL
  inf <- sort(inf,by=~seqnames + start + end)
}

if(grepl(pattern = '[.]bed',ref_file)){
  ref <- unlist(blocks(ref))
  ref$gene_id <- gsub(';.*','',names(ref))
  ref$transcript_id <- gsub('.*;','',names(ref))
  ref$type <- 'exon'
  names(ref) <- NULL
  ref <- sort(ref,by=~seqnames + start + end)
}

ref <- cleanMcols(ref)
inf <- cleanMcols(inf)

ref <- ref[ref$type=='exon',]
inf <- inf[inf$type=='exon',]

## filter gene with multiple strands in inf RTD
strand_idx <- elementNROWS(range(GenomicRanges::split(inf,inf$gene_id)))
strand_idx <- names(strand_idx)[strand_idx>1]
inf <- inf[!(inf$gene_id %in% strand_idx),]

## add prefix to gene ids
ref$gene_id <- paste0(prefix_ref,'_',ref$gene_id)
ref$transcript_id <- paste0(prefix_ref,'_',ref$transcript_id)
ref$source <- prefix_ref

inf$gene_id <- paste0(prefix_inf,'_',inf$gene_id)
inf$transcript_id <- paste0(prefix_inf,'_',inf$transcript_id)
inf$source <- prefix_inf

message('Step 2: Identify novel transcript in the inference RTD')
inf_trans_list <- list()

message('  -> Process multi-exon transcripts')
####################################################
###---> Multi-exon transcripts

## generate unique labels of SJ chain
intron_ref <- exon2intron(ref)
intron_ref <- sort(intron_ref,by = ~ seqnames + start + end)
intron_ref_sj <- intron_ref$label
intron_ref_trans <- base::split(intron_ref_sj,intron_ref$transcript_id)
intron_ref_trans <- sapply(intron_ref_trans, function(x) paste(x,collapse = ';'))

intron_inf <- exon2intron(inf)
intron_inf <- sort(intron_inf,by = ~ seqnames + start + end)
intron_inf_sj <- intron_inf$label
intron_inf_trans <- base::split(intron_inf_sj,intron_inf$transcript_id)
intron_inf_trans <- sapply(intron_inf_trans, function(x) paste(x,collapse = ';'))

## add source the ref rtd
trans_add_s <- names(intron_ref_trans)[intron_ref_trans %in% intron_inf_trans]
ref$source[ref$transcript_id %in% trans_add_s] <- paste0('SJ:',prefix_ref,'&',prefix_inf)

## get the transcript names with novel SJ labels
idx <- which(!(intron_inf_sj %in% intron_ref_sj))
trans_multiexon_novel <- unique(intron_inf[idx]$transcript_id)

## multi-exon transcripts in inf which are not novel
trans_multiexon_filtered <- setdiff(names(intron_inf_sj),trans_multiexon_novel)
inf_trans_list <- c(inf_trans_list,trans_multiexon_filtered=list(trans_multiexon_filtered))
inf_trans_list <- c(inf_trans_list,trans_multiexon_novel=list(trans_multiexon_novel))

message('  -> Process mono-exon transcripts')
## get mono exon transcripts
ref_mono <- getMonoExonTrans(ref)
inf_mono <- getMonoExonTrans(inf)

## check whether inf mono exon trans has overlap of exons of ref, introic is excluded.
trans_mono <- unique(inf_mono$transcript_id)
hits <- suppressWarnings(findOverlaps(query = inf_mono,
                                      subject = ref,ignore.strand=FALSE))

## get novel mono exon trans with no overlap
trans_overlap <- inf_mono[queryHits(hits),]$transcript_id
trans_monoexon_novel <- setdiff(trans_mono,trans_overlap)
trans_monoexon_filtered <- setdiff(trans_mono,trans_monoexon_novel)

inf_trans_list <- c(inf_trans_list,trans_monoexon_novel=list(trans_monoexon_novel))
inf_trans_list <- c(inf_trans_list,trans_monoexon_filtered=list(trans_monoexon_filtered))

## novel transcripts in inf rtd, multi-exon + mono-exon
trans_novel<- union(trans_multiexon_novel,trans_monoexon_novel)
inf2keep <- inf[inf$transcript_id %in% trans_novel]

message('Step 3: Merge inf to ref RTD and rename transcript and gene ids')

rtd <- suppressWarnings(c(ref,inf2keep))
rtd <- sort(rtd,by=~seqnames+start+end)
rtd$gene_id_old <- rtd$gene_id
rtd$transcript_id_old <- rtd$transcript_id

if(check_canonical){
  if(is.null(genome_file)){
    message('Genome is not provide, Canonical SJ is not checked')
  } else {
    genome <- import(genome_file)
    motif <- getIntronMotif(gr = rtd, genome = genome)
    trans <- unique(motif$transcript_id[motif$canonical=='No'])
    if(length(trans) > 0)
      rtd <- rtd[!(rtd$transcript_id %in% trans)] else rtd <- rtd
  }
}


### Gene initial rename, overlap between transcripts and gene collapsed region
gene_overlap <- getOverlapRange(rtd)
trans_range <- getTransRange(rtd)
hits <- findOverlaps(query = trans_range,subject = gene_overlap)

### Get the initial gene-transcript mapping table
mapping <- data.frame(seqnames=seqnames(trans_range[queryHits(hits)]),
                      transcript_id_old=trans_range[queryHits(hits)]$transcript_id,
                      gene_id_group=gene_overlap[subjectHits(hits)]$gene_id)
mapping <- mapping[!duplicated(mapping),]
mapping <- mapping[order(mapping$gene_id_group,decreasing = F),]
mapping$gene_i <- rep(1:length(unique(mapping$gene_id_group)),
                      rle(mapping$gene_id_group)$lengths)
pad_n <- length(unique(mapping$gene_i))
pad_n <- nchar(format(pad_n,scientific = FALSE))
mapping$gene_id <- paste0(mapping$seqnames,'G',
                          stringi::stri_pad(mapping$gene_i,width = pad_n,pad = '0'))
mapping$transcript_i <- sequence(rle(mapping$gene_id)$lengths)
mapping$transcript_id <- paste0(mapping$gene_id,'.',mapping$transcript_i)
rownames(mapping) <- mapping$transcript_id_old

## rename
rtd$transcript_id <- mapping[rtd$transcript_id_old,'transcript_id']
rtd$gene_id <- mapping[rtd$transcript_id_old,'gene_id']

##### check intronic
# overlap between trans_range and intron disjoin
message('  -> Process intronic genes')
trans_range <- getTransRange(rtd)
introns <- exon2intron(exons = rtd)
introns_s <- GenomicRanges::split(introns,introns$gene_id)
introns_d <- unlist(disjoin(introns_s))
introns_d$gene_id <- names(introns_d)

hits <- findOverlaps(trans_range,introns_d,type = 'within')
if(NROW(hits)>0){
  intronic <- trans_range[queryHits(hits)]
  gene_overlap_intronic <- getOverlapRange(intronic)
  trans_range_intronic <- getTransRange(intronic)
  hits <- findOverlaps(query = trans_range_intronic,subject = gene_overlap_intronic)
  mapping_intronic  <- data.frame(
    seqnames=seqnames(trans_range_intronic[queryHits(hits)]),
    transcript_id_old=trans_range_intronic[queryHits(hits)]$transcript_id,
    gene_id_group=gene_overlap_intronic[subjectHits(hits)]$gene_id)
  
  mapping_intronic <- mapping_intronic[order(mapping_intronic$gene_id_group,
                                             decreasing = F),]
  mapping_intronic$gene_i <- rep(1:length(unique(mapping_intronic$gene_id_group)),
                                 rle(mapping_intronic$gene_id_group)$lengths)
  n <- max(mapping_intronic$gene_i)
  mapping_intronic$gene_i <- mapping_intronic$gene_i/10^ceiling(log10(n)+1)
  rownames(mapping_intronic) <- mapping_intronic$transcript_id_old
  
  rownames(mapping) <- mapping$transcript_id
  idx <- mapping_intronic$transcript_id_old
  mapping[idx,'gene_i'] <- mapping[idx,'gene_i']+mapping_intronic[idx,'gene_i']
  mapping$gene_i <- as.integer(factor(mapping$gene_i))
  
  mapping$gene_id <- paste0(mapping$seqnames,'G',
                            stringi::stri_pad(mapping$gene_i,width = pad_n,pad = '0'))
  mapping <- mapping[order(mapping$gene_id,decreasing = F),]
  mapping$transcript_i <- sequence(rle(mapping$gene_id)$lengths)
  mapping$transcript_id <- paste0(mapping$gene_id,'.',mapping$transcript_i)
  rownames(mapping) <- mapping$transcript_id_old
  rtd$transcript_id <- mapping[rtd$transcript_id_old,'transcript_id']
  rtd$gene_id <- mapping[rtd$transcript_id_old,'gene_id']
  
  trans_range <- getTransRange(rtd)
}
# trans_range$gene_id <- gsub('[.].*','',trans_range$transcript_id)
trans_range$gene_id <- sub("\\.[^.]*$", "", trans_range$transcript_id)


message('  -> Process chimeric check')
###---> find overlap between transcript ranges and transcript ranges
hits <- findOverlaps(query = trans_range,subject = trans_range)
idx <- which(queryHits(hits)!=subjectHits(hits))
hits <- hits[idx,]

## get the overlaps < chimeric_tolerance
gr1 <- trans_range[queryHits(hits)]
gr2 <- trans_range[subjectHits(hits)]
overlaps <- pintersect(gr1,gr2)
p1 <- width(overlaps)/width(gr1)
p2 <- width(overlaps)/width(gr2)
idx <- which(p1 < chimeric_tolerance & p2 < chimeric_tolerance)
overlaps_pass <- overlaps[idx]
overlaps_pass$from <- gr1$transcript_id[idx]
overlaps_pass$to <- gr2$transcript_id[idx]
overlaps_pass$transcript_id <- NULL
overlaps_pass$hit <- NULL
overlaps_pass <- unique(overlaps_pass)

## do second check, if the overlaps inside in other transcritps
trans <- union(gr1$transcript_id[idx],gr2$transcript_id[idx])
genes <- unique(sub('\\..[^\\.]*$', '', trans))
trans_range_pass <- trans_range[trans_range$gene_id %in% genes]
trans_range_pass_remain <-
  trans_range_pass[!(trans_range_pass$transcript_id %in% trans)]
hits_pass <- findOverlaps(query = overlaps_pass,trans_range_pass_remain,
                          type = 'within')

idx <- unique(queryHits(hits_pass))
overlaps_pass <- overlaps_pass[-idx,]

overlaps_pass_split <- GenomicRanges::split(overlaps_pass,overlaps_pass$gene_id)
## if multiple overlap, get range and check range length agaist transcript length
if(any(elementNROWS(overlaps_pass_split)>1)){
  overlaps_pass_split <- unlist(range(overlaps_pass_split))
  overlaps_pass_split$gene_id <- names(overlaps_pass_split)
  # names(overlaps_pass_split) <- NULL
  overlaps_pass_split <- overlaps_pass_split[overlaps_pass$gene_id]
  overlaps_pass_split$from <- overlaps_pass$from
  overlaps_pass_split$to <- overlaps_pass$to
  
  #check again if the region is more than > tolorance
  trans_l <- width(trans_range)
  names(trans_l) <- trans_range$transcript_id
  overlaps_pass_split$from_overlap <- width(overlaps_pass_split)/trans_l[overlaps_pass_split$from]
  overlaps_pass_split$to_overlap <- width(overlaps_pass_split)/trans_l[overlaps_pass_split$to]
  
  idx_filter <- which(overlaps_pass_split$to_overlap >= chimeric_tolerance |
                        overlaps_pass_split$from_overlap >= chimeric_tolerance)
  if(length(idx_filter)>0){
    genes2filer <- unique(overlaps_pass_split$gene_id[idx_filter])
    overlaps_pass <- overlaps_pass[!(overlaps_pass$gene_id %in% genes2filer)]
  } else {
    overlaps_pass <- overlaps_pass
  }
}

trans_range_pass <- trans_range_pass[trans_range_pass$gene_id %in% overlaps_pass$gene_id]

## repeat the overlap to match the overlap to all the transcripts in the gene
overlaps_pass <- GenomicRanges::split(overlaps_pass,overlaps_pass$gene_id)
overlaps_pass <- reduce(overlaps_pass)
overlaps_pass <- unlist(overlaps_pass)
overlaps_pass <- overlaps_pass[trans_range_pass$gene_id]
overlaps_pass$gene_id <- trans_range_pass$gene_id
overlaps_pass$transcript_id <- trans_range_pass$transcript_id

## check agian if any gene in the middle of the transcript ranges
idx2rm <- which(start(overlaps_pass) > start(trans_range_pass) &
                  end(overlaps_pass) < end(trans_range_pass))
if(length(idx2rm) > 0){
  genes2rm <- unique(trans_range_pass$gene_id[idx2rm])
  trans_range_pass <- trans_range_pass[!(trans_range_pass$gene_id %in% genes2rm)]
  overlaps_pass <- overlaps_pass[trans_range_pass$gene_id]
}

if(NROW(overlaps_pass) > 0 ){
  ## take away overlap from all the transcripts
  trans_range_pass_cut <- psetdiff(trans_range_pass,overlaps_pass)
  mcols(trans_range_pass_cut) <- mcols(trans_range_pass)
  
  gene_overlap_chemric <- getOverlapRange(trans_range_pass_cut)
  hits <- findOverlaps(trans_range_pass_cut,gene_overlap_chemric)
} else {
  hits <- NULL
}

if(NROW(hits)>0){
  mapping_chemric  <- data.frame(
    seqnames=seqnames(trans_range_pass_cut[queryHits(hits)]),
    transcript_id_old=trans_range_pass_cut[queryHits(hits)]$transcript_id,
    gene_id_group=gene_overlap_chemric[subjectHits(hits)]$gene_id)
  
  mapping_chemric <- mapping_chemric[order(mapping_chemric$gene_id_group,
                                           decreasing = F),]
  mapping_chemric$gene_i <- rep(1:length(unique(mapping_chemric$gene_id_group)),
                                rle(mapping_chemric$gene_id_group)$lengths)
  n <- max(mapping_chemric$gene_i)
  mapping_chemric$gene_i <- mapping_chemric$gene_i/10^ceiling(log10(n)+1)
  rownames(mapping_chemric) <- mapping_chemric$transcript_id_old
  
  rownames(mapping) <- mapping$transcript_id
  idx <- mapping_chemric$transcript_id_old
  mapping[idx,'gene_i'] <- mapping[idx,'gene_i']+mapping_chemric[idx,'gene_i']
  mapping$gene_i <- as.integer(factor(mapping$gene_i))
  mapping$gene_id <- paste0(mapping$seqnames,'G',
                            stringi::stri_pad(mapping$gene_i,width = pad_n,pad = '0'))
  mapping$transcript_i <- sequence(rle(mapping$gene_id)$lengths)
  mapping$transcript_id <- paste0(mapping$gene_id,'.',mapping$transcript_i)
  rownames(mapping) <- mapping$transcript_id_old
  rtd$transcript_id <- mapping[rtd$transcript_id_old,'transcript_id']
  rtd$gene_id <- mapping[rtd$transcript_id_old,'gene_id']
}

rtd$gene_id_old <- NULL
rtd$transcript_id_old <- NULL
if(!is.null(prefix_gene)){
  rtd$gene_id <- paste0(prefix_gene,'_',rtd$gene_id)
  rtd$transcript_id <- paste0(prefix_gene,'_',rtd$transcript_id)
}

rtd <- sort(rtd,by = ~ seqnames + start + end)

message('Step 4: Summary the merged RTD')

###---> transcript per gene
n <- as.character(c(1:20,'>20'))
diversity <- lapply(list(ref=ref,inf=inf,merged=rtd), function(gr){
  trans2gene <- data.frame(TXNAME=gr$transcript_id,GENEID=gr$gene_id)
  trans2gene <- trans2gene[!duplicated(trans2gene),]
  idx <- table(table(trans2gene$GENEID))
  if(max(as.numeric(names(idx))) > 20){
    idx1 <- idx[as.numeric(names(idx))<=20]
    idx2 <- sum(idx[as.numeric(names(idx))>20])
    names(idx2) <- '>20'
    idx <- c(idx1,idx2)
  }
  num <- rep(0,length(n))
  names(num) <- n
  num[names(idx)] <- idx
  num
})

transpergene <- do.call(cbind,diversity)
transpergene <- data.frame(TransPerGene=rownames(transpergene),transpergene,row.names = NULL)
colnames(transpergene) <- c('TransPerGene',prefix_ref,prefix_inf,'merged')
write.csv(transpergene,
          file=file.path(data_dir,paste0('htsm_rtd_',prefix_ref,'_', prefix_inf,
                                         '_transcript_per_gene_number.csv')),
          row.names = F)

data2plot <- transpergene[,-1]
rownames(data2plot) <- transpergene[,1]
data2plot <- t(data2plot)

png(file.path(data_dir,
              paste0('htsm_rtd_',prefix_ref,'_',prefix_inf,'_transcript_per_gene_number.png')),
    width = 10,height = 5,res = 300,units = 'in')
xx <- barplot(data2plot,beside = TRUE,
              col =  c("#999999", "#E69F00", "#56B4E9"),
              ylim=c(0,1.2*max(data2plot)),
              xlab = 'Transcripts in a gene',
              ylab = 'Genen number',
              main='Transcript number per gene',
              space=c(0,0.5))
## Add text at top of bars
text(x = xx, y = data2plot, label = data2plot, cex = 0.7, col = "red",
     srt=-90,adj=c(1.1,0.5))
legend("topright",
       legend = rownames(data2plot),
       fill = c("#999999", "#E69F00", "#56B4E9"))
dev.off()


pdf(file.path(data_dir,
              paste0('htsm_rtd_',prefix_ref,'_',prefix_inf,'_transcript_per_gene_number.pdf')),
    width = 10,height = 5)
xx <- barplot(data2plot,beside = TRUE,
              col =  c("#999999", "#E69F00", "#56B4E9"),
              ylim=c(0,1.2*max(data2plot)),
              xlab = 'Transcripts in a gene',
              ylab = 'Genen number',
              main='Transcript number per gene',
              space=c(0,0.5))

## Add text at top of bars
text(x = xx, y = data2plot, label = data2plot, cex = 0.7, col = "red",
     srt=-90,adj=c(1.1,0.5))
legend("topright",
       legend = rownames(data2plot),
       fill = c("#999999", "#E69F00", "#56B4E9"))
dev.off()

###---> basic statistics
stat_merged <- rtdSummary(gr = rtd)
stat_ref <- rtdSummary(gr = ref)
stat_inf <- rtdSummary(gr = inf)
stat <- cbind(stat_ref,stat_inf,stat_merged)
colnames(stat) <- c(prefix_ref,prefix_inf,'merged')

write.csv(stat,file=file.path(data_dir,paste0('htsm_rtd_',prefix_ref,'_', prefix_inf,'_merged_summary.csv')))
message('Step 5: Save the results')

export(object = rtd,
       con = file.path(data_dir,paste0('htsm_rtd_',prefix_ref,'_', prefix_inf,'_merged.gtf')))
# export_gtf(gr = rtd,file2save = file.path(data_dir,paste0('htsm_rtd_',prefix_ref,'_', prefix_inf,'_merged.gtf')))
# save(rtd,file = file.path(data_dir,paste0('htsm_rtd_',prefix_ref,'_', prefix_inf,'_merged.RData')))
# save(inf_trans_list,file=file.path(data_dir,paste0(prefix_inf,'_trans_summary.RData')))

#######################
###---> save to bed file
exon <- rtd
mapping <- data.frame(TXNAME = exon$transcript_id,GENEID = exon$gene_id)
mapping <- mapping[!duplicated(mapping),]
rownames(mapping) <- mapping$TXNAME

exon <- GenomicRanges::split(exon,exon$transcript_id)
bed <- rtracklayer::asBED(exon)

blocks <- bed$blocks
bed$blocks <- NULL
strand_info <- as.character(strand(bed))
strand_info[strand_info == "*"] <- NA

df <- data.frame(seqnames = seqnames(bed),
                 start = start(bed) - 1,
                 end = end(bed),
                 names = paste0(mapping[bed$name,'GENEID'],';',bed$name),
                 score = 0,
                 strand = strand_info,
                 thickStart = start(bed) - 1,
                 thickEnd = end(bed),
                 itemRgb = 1,
                 blockCount = elementNROWS(blocks),
                 blockSizes = unlist(lapply(width(blocks), paste, collapse = ","), use.names = FALSE),
                 blockStarts = unlist(lapply(start(blocks) - 1, paste, collapse = ","), use.names = FALSE))
rownames(df) <- NULL

message('Export bed file ...')
write.table(df, file = file.path(data_dir,paste0('htsm_rtd_',prefix_ref,'_', prefix_inf,'_merged.bed')),
            sep = "\t", col.names = FALSE,
            row.names = FALSE, quote = FALSE, na = ".")

write.csv(mapping, file = file.path(data_dir,paste0('htsm_rtd_',prefix_ref,'_', prefix_inf,'_merged_trans_gene_mapping.csv')),
          row.names = FALSE)

message('Done: ',Sys.time())
end.time <- Sys.time()
time.taken <- end.time - start.time
message(paste('Time for analysis:',round(time.taken,3),attributes(time.taken)$units),'\n')
```

<a href='#table-of-contents'>Go back to Table of contents</a>

Barley PanBaRT20
-----------


### Construction of linear pan-genome

The PSVCP method (https://github.com/wjian8/psvcp_v1.01) was used to generate linear barley pan-gemoe. Precisely, the Morex reference genome served as the reference, and through an iterative process, 19 additional genomes were progressively integrated to build a comprehensive linear pan-genome. The order of the integration: 

| Order | Barley   reference genome                                                         |
|-------|-----------------------------------------------------------------------------------|
| 0     | 210316_Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.fasta.gz               |
| 1     | 220812_FT11_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz                 |
| 2     | 220214_Golden_Promise_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz       |
| 3     | 220503_HOR_8148_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz             |
| 4     | 220816_RGT_Planet_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz           |
| 5     | 220613_HOR_21599_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz            |
| 6     | 220518_HOR_13821_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz            |
| 7     | 220613_OUN333_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz               |
| 8     | 220519_HOR_3365_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz             |
| 9     | 220503_HOR_9043_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz             |
| 10    | 220613_HOR_10350_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz            |
| 11    | 220519_HOR_13942_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz            |
| 12    | 220816_Akashinirki_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz          |
| 13    | 220503_HOR_7552_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz             |
| 14    | 220816_Du_Li_Huang_ZDM01467_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz |
| 15    | 220411_HOR_3081_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz             |
| 16    | 211124_Barke_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz                |
| 17    | 220810_Igri_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz                 |
| 18    | 220809_Hocket_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz               |
| 19    | 220809_Chiba_ZDM02064_pseudomolecules_and_unplaced_contigs_CPclean.fasta.gz       |


```
python3 $path_of_the_pipeline/Construct_pan_and_Call_sv.py 
    genome_gff_dir 
    genome_list 
    -fqd fq_dir 
    -o population_hmp

```


### Map GsRTDs to linear pan-genome

Minimap2 v2.24 (Li, 2018) was used for the gene/transcript mapping: 

  -   Step 1: Mapping of GsRTD gene sequences (encompassing exons and introns) to the linear pan-genome.
  -   Step 2: Mapping of GsRTD transcript sequences to the specific region on the linear pan-genome where their corresponding genes are located. This two-step process ensures that transcripts from the same gene are not mapped to different gene loci, preventing potential errors. 

```
# genome_fasta: the linear pan-genome fasta file
# trans_fasta: the gene or transcript sequence fasta files in above two steps.
# output_sam: the output sam file. 

minimap2 \
-t 16 \
-G 15,000 \
-L \
--secondary=no \
--MD \
-ax splice:hq -uf $genome_fasta $trans_fasta > $output_sam

prefix=${output_sam%.*}
samtools view -Sb -F 2048 ${prefix}.sam > ${prefix}.bam
samtools sort ${prefix}.bam > ${prefix}_sorted.bam

bedtools bamtobed -bed12 -i ${prefix}_sorted.bam > ${prefix}.bed
bedToGenePred ${prefix}.bed ${prefix}.genepred
genePredToGtf "file" ${prefix}.genepred ${prefix}.gtf
```

### Gene and trasncript association between PanBaRT20 and GsRTDs

In-house R scripts were used to generate the association and the look-up table between PanBaRT20 and GsRTDs. 

  - Step 1: Based on the GsRTD transcript mapping on the linear pan-genome, the transcripts were grouped by overlap-layout-consensus. (1) The overlapped transcripts were assigned with the same gene ID. (2) The multiple-exon transcripts with identical intron combinations but different TSS or TES were merged to be one PanBaRT20 transcript and the longest TSS and TES were used as the starting and ending point of this transcript. The overlapped mono-exon transcripts were all merged into one PanBaRT20 transcript. 
  - Step 2: Refine gene IDs in the PanBaRT20. (1) If two overlapped transcript sets only overlap 5% of their lengths, they were treated as different genes and assigned with different IDs. (2) If an overlapped transcript set was included entriely within the intron region of another transcript set, they were treated as different genes and assigned with different IDs. 

Thus, each transcripts in PanBaRT20 represent a set of GsRTD transcripts and we used this information to generate a look-up table between PanBaRT20 and GsRTDs. 

```
#######################################
##---> R code of step 1

# gtf_file: the output from previous minimap2 result, which mapped the GsRTD transcripts to the linear pan-genome.


library(rtracklayer)
gr <- import(gtf_file)
###############
###---> mutli exon transcripts
multi_trans <- getMultiExonTrans(gr)
intron_chain <- getIntronChain(multi_trans)
multi_trans$group <- intron_chain[multi_trans$transcript_id] 
multi_trans_group <- GenomicRanges::split(multi_trans$transcript_id,multi_trans$group)
multi_trans_group <- lapply(multi_trans_group, unique)

multi_trans_split <- GenomicRanges::split(multi_trans,multi_trans$group)
multi_trans_split <- reduce(multi_trans_split)
multi_trans <- unlist(multi_trans_split)
multi_trans$type <- 'exon'
multi_trans$gene_id <- names(multi_trans)
multi_trans$transcript_id <- names(multi_trans)
names(multi_trans) <- NULL  
# export(multi_trans,'multi_trans.gtf')

###---> mono exon transcripts
mono_trans <- getMonoExonTrans(gr)
mcols(mono_trans) <- mcols(mono_trans)[,c('type','gene_id','transcript_id')]
mono_loci <- getOverlapRange(mono_trans)
hits <- findOverlaps(query = mono_trans,subject = mono_loci,type = 'within')
mono_trans <- mono_trans[hits@from]
mono_trans$group <- mono_loci[hits@to]$gene_id
mono_trans <- unique(mono_trans)

mono_trans_group <- split(mono_trans$transcript_id,mono_trans$group)
mono_trans_group <- lapply(mono_trans_group, unique)

mono_trans_split <- GenomicRanges::split(mono_trans,mono_trans$group)
mono_trans_split <- reduce(mono_trans_split)
mono_trans <- unlist(mono_trans_split)
mono_trans$type <- 'exon'
mono_trans$gene_id <- names(mono_trans)
mono_trans$transcript_id <- names(mono_trans)
names(mono_trans) <- NULL 

trans_group <- c(mono_trans_group,multi_trans_group)
trans_group <- sapply(trans_group, function(x) paste0(unique(x),collapse = ';'))
names(trans_group) <- paste0('trans-',names(trans_group))

gr_collapsed <- c(multi_trans,mono_trans)
gr_collapsed$transcript_id <- paste0('trans-',gr_collapsed$transcript_id)
gr_collapsed <- sort(gr_collapsed,by = ~ seqnames + start + end)
gr_collapsed$source_trans <- trans_group[gr_collapsed$transcript_id]

#######################################
##---> R code of step 2

gr <- gr_collapsed
## record original gene and transcript id
gr$gene_id0 <- gr$gene_id
gr$transcript_id0 <- gr$transcript_id
gr$observation <- NA

trans_range <- getTransRange(gr)
hits <- findOverlaps(query = trans_range,subject = trans_range)
idx <- which(queryHits(hits)!=subjectHits(hits))
hits <- hits[idx,]

## get the overlaps < chimeric_tolerance
message('Refine chimeric genes')
gr1 <- trans_range[queryHits(hits)]
gr2 <- trans_range[subjectHits(hits)]
overlaps <- pintersect(gr1,gr2)
overlaps$pair1 <- gr1$transcript_id
overlaps$pair2 <- gr2$transcript_id

p1 <- width(overlaps)/width(gr1)
p2 <- width(overlaps)/width(gr2)
idx <- which(p1 < chimeric_tolerance & p2 < chimeric_tolerance)

if(length(idx) > 0){
  chemric <- overlaps[idx]
  
  # names(chemric) <- chemric$transcript_id
  pair1 <- GenomicRanges::split(chemric,chemric$pair1)
  pair1 <- reduce(pair1)
  pair1 <- unlist(pair1)
  pair1$transcript_id <- names(pair1)
  names(pair1) <- NULL
  
  #chr1H:99,372,722-99,380,409
  pair2 <- GenomicRanges::split(chemric,chemric$pair2)
  pair2 <- reduce(pair2)
  pair2 <- unlist(pair2)
  pair2$transcript_id <- names(pair2)
  names(pair2) <- NULL
  
  chemric <- c(pair1,pair2)
  chemric <- sort(chemric,by = ~ seqnames + start + end)
  gr_chemric <- gr[gr$transcript_id %in% chemric$transcript_id]
  gr_no_chemric <- gr[!(gr$transcript_id %in% chemric$transcript_id)]
  # export(chemric,'chemricchemric.gtf')
  
  chemric <- GenomicRanges::split(chemric,chemric$transcript_id)
  chemric <- chemric[gr_chemric$transcript_id]
  gr_chemric_cut <- psetdiff(gr_chemric,chemric)
  gr_chemric_cut <- unlist(gr_chemric_cut)
  gr_chemric_cut$transcript_id <- names(gr_chemric_cut)
  names(gr_chemric_cut) <- NULL
  gr_chemric_cut$gene_id <- gr_chemric_cut$transcript_id
  gr_chemric_cut$type <- 'exon'
  
  gr_cut <- c(gr_chemric_cut,gr_no_chemric)
  gr_cut <- sort(gr_cut,by = ~ seqnames + start + end)
  gr_cut <- assignGeneID(gr = gr_cut)
  mapping <- data.frame(transcript_id=gr_cut$transcript_id,
                        gene_id=gr_cut$gene_id)
  mapping <- unique(mapping)
  rownames(mapping) <- mapping$transcript_id
  gr$gene_id <- mapping[gr$transcript_id,'gene_id']
}

##################
###---> intronic
message('Refine intronic genes')
gr_intronic <- getIntronicGene(query = gr,subject = gr)
if(NROW(gr_intronic)>0){
  gr_intronic <- assignGeneID(gr = gr_intronic)
  gr_intronic$gene_id <- paste0('intronic-',gr_intronic$gene_id)
  mapping <- data.frame(transcript_id=gr_intronic$transcript_id,
                        gene_id=gr_intronic$gene_id)
  mapping <- unique(mapping)
  rownames(mapping) <- mapping$transcript_id
  
  idx <- which(gr$transcript_id %in% gr_intronic$transcript_id)
  gr$observation[idx] <- 'intronic'
  gr$gene_id[idx] <- mapping[gr$transcript_id[idx],'gene_id']
}

#########
##---> generate new gene ids
message('Generate new gene ids')
gr <- sort(gr,by = ~ seqnames + start + end)
genes <- unique(gr$gene_id)

pad_n <- length(genes)
pad_n <- nchar(format(pad_n,scientific = FALSE))

genes_n <- stringr::str_pad(string = 1:length(genes),width = pad_n,pad = '0')
names(genes_n) <- genes

gene_id <- paste0(prefix,seqnames(gr),genes_n[gr$gene_id])
gr$gene_id <- gene_id

mapping <- data.frame(transcript_id=gr$transcript_id,
                      gene_id=gene_id)
mapping <- unique(mapping)
mapping <- mapping[order(mapping$gene_id),]
n <- sequence(rle(mapping$gene_id)$lengths)
transcript_id <- paste0(mapping$gene_id,'.',n)
names(transcript_id) <- mapping$transcript_id
gr$transcript_id <- transcript_id[gr$transcript_id]

#############################
###---> identify multi-exon and mono-exon transcripts
gr_multi <- getMultiExonTrans(gr)
idx <- which(gr$transcript_id %in% gr_multi$transcript_id)
if(length(idx) > 0)
  gr$observation[idx] <- paste0(gr$observation[idx],';multi-exon')

gr_mono <- getMonoExonTrans(gr)
idx <- which(gr$transcript_id %in% gr_mono$transcript_id)
if(length(idx) > 0)
  gr$observation[idx] <- paste0(gr$observation[idx],';mono-exon')

###---> identify exonic transcripts
trans_range <- getTransRange(gr_mono)
hits <- findOverlaps(trans_range,gr,type = 'within')
mapping <- data.frame(from=trans_range$transcript_id[hits@from],
                      to=gr$transcript_id[hits@to])
mapping <- unique(mapping)
trans <- unique(mapping$from[mapping$from != mapping$to])
idx <- which(gr$transcript_id %in% trans)
if(length(idx) > 0)
  gr$observation[idx] <- paste0(gr$observation[idx],';exonic-trans')

gr$observation <- gsub('NA;','',gr$observation)

export(gr,'PanBaRT20.gtf')
write.table(mapping2, file = "GsRTD_and_PanBaRT20_match.tsv", row.names=FALSE, sep="\t")
```


Expression analysis
----------

### Transcript quantification

Salmon v1.10.1 was used for transcript quantification.

#### Salmon index

```
# fasta_file: transcript sequence file
# index_dir: output index direcotry

salmon index -t ${fasta_file} \
-i ${index_dir} \
-p 12 \
-k 31 \
--keepDuplicates
```

#### Salmon quantification

```
# index_dir: output index direcotry
# read1: fastq file of paired-end read1
# read2: fastq file of paired-end read2
# output_folder: output directory

salmon quant \
-i $index_dir \
-l ISR \
-1 $read1 \
-2 $read2 \
-p 12 \
-o $output_folder \
--seqBias \
--posBias \
--gcBias \
--validateMappings
```

### Differential expression analysis

The 3D RNA-seq App was used to perform differential gene expression, differential gene alternative splicing, differential transcript expression and differential transcript usage analysis. 

References
----------
  -   Chen,J., Tang,X., Ren,C., Wei,B., Wu,Y., Wu,Q., and Pei,J. (2018) Full-length transcriptome sequences and the identification of putative genes for flavonoid biosynthesis in safflower. BMC Genomics, 19, 1–13.
  -   Coulter,M., Entizne,J.C., Guo,W., Bayer,M., Wonneberger,R., Milne,L., Schreiber,M., Haaning,A., Muehlbauer,G., McCallum,N., Fuller,J., Simpson,C., Stein,N., Brown,J.W.S., Waugh,R., and   -     -   Zhang,R. (2022) BaRTv2: A highly resolved barley reference transcriptome for accurate transcript‐specific RNA‐seq quantification. Plant J., 111, 1183–1202.
  -   Dobin,A., Davis,C.A., Schlesinger,F., Drenkow,J., Zaleski,C., Jha,S., Batut,P., Chaisson,M., and Gingeras,T.R. (2013) STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29, 15–21.
  -   Dobin,A. and Gingeras,T.R. (2015) Mapping RNA-seq Reads with STAR. Curr. Protoc. Bioinforma., 51, 11.14.1-11.14.19.
Entizne,J.C., Guo,W., Calixto,C.P., Spensley,M., Tzioutziou,N.,   -   Zhang,R., and Brown,J.W. (2020) TranSuite: a software suite for accurate translation and characterization of transcripts. bioRxiv, 2020.12.15.422989.
  -   Guo,W., Tzioutziou,N.A., Stephen,G., Milne,I., Calixto,C.P.G., Waugh,R., Brown,J.W.S., and Zhang,R. (2021) 3D RNA-seq: a powerful and flexible tool for rapid and accurate differential expression and alternative splicing analysis of RNA-seq data for biologists. RNA Biol., 18, 1574–1587.
  -   Jayakodi,M., Padmarasu,S., Haberer,G., Bonthala,V.S., Gundlach,H., Monat,C., Lux,T., Kamal,N., Lang,D., Himmelbach,A., Ens,J., Zhang,X.Q., Angessa,T.T., Zhou,G., Tan,C., Hill,C., Wang,P.,   -   Schreiber,M., Boston,L.B., et al. (2020) The barley pan-genome reveals the hidden legacy of mutation breeding. Nature, 588, 284–289.
  -   Kuo,R.I., Cheng,Y., Zhang,R., Brown,J.W.S., Smith,J., Archibald,A.L., and Burt,D.W. (2020) Illuminating the dark side of the human transcriptome with long read transcript sequencing. BMC Genomics, 21, 1–22.
  -   Li,H. (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34, 3094–3100.
Patro,R., Duggal,G., Love,M.I., Irizarry,R.A., and Kingsford,C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417–419.
  -   Pertea,M., Pertea,G.M., Antonescu,C.M., Chang,T.-C., Mendell,J.T., and Salzberg,S.L. (2015) StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat. Biotechnol., 33, 290–295.
  -   Quevillon,E., Silventoinen,V., Pillai,S., Harte,N., Mulder,N., Apweiler,R., and Lopez,R. (2005) InterProScan: Protein domains identifier. Nucleic Acids Res., 33, W116.
  -   Shao,M. and Kingsford,C. (2017) Accurate assembly of transcripts through phase-preserving graph decomposition. Nat. Biotechnol., 35, 1167–1169.
  -   Wang,J., Yang,W., Zhang,S., Hu,H., Yuan,Y., Dong,J., Chen,L., Ma,Y., Yang,T., Zhou,L., Chen,J., Liu,B., Li,C., Edwards,D., and Zhao,J. (2023) A pangenome analysis pipeline provides insights into functional gene identification in rice. Genome Biol., 24, 1–22.


<a href='#table-of-contents'>Go back to Table of contents</a>

</div>
