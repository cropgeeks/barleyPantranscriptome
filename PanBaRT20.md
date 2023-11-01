---
title: "README"
author: "Wenbin Guo"
date: "2023-11-01"
output: html_document
---


# PanBaRT20 analysis pipeline

<hr>

Table of contents
-----------------

-   [Description](#description)
-   [Short-read RTD construction](#short-read-rtd-construction)
-   [Long-read RTD construction](#long-read-rtd-construction)
-   [Genotype-specific RTD construction](#genotype-specific-rtd-construction)
-   [Barley Pan-RTD construction](#barley-pan-rtd-construction)
-   [References](#references)

<div align="justify">



Description
-----------


Short-read RTD construction
-----------

### Data pre-processing

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

### Data pre-processing
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

### Read mapping

FLNC reads were mapped to reference genome by using Minimap2 v2.24 (Li, 2018). The maximum intron size was set to 15,000. 

```

# genome: reference genome fasta file
# input: FLNC read
# output: output directory of the mapping

minimap2 -ax splice:hq -uf -G $maxintron ${genome} ${input} -o ${output}

```

### Transcript assembly

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


### Splice junction analysis


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


### Transcript start and end analysis


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


### Filter transcripts with low-quality SJs or TSS/TES 

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



Genotype-specific RTD construction
----------- 


Barley Pan-RTD construction
----------- 


References
----------

-   Bray,N.L., Pimentel,H., Melsted,P., and Pachter,L. (2016) Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol., 34, 525–527.
-   Calixto,C.P.G., Guo,W., James,A.B., Tzioutziou,N.A., Entizne,J.C., Panter,P.E., Knight,H., Nimmo,H.G., Zhang,R., and Brown,J.W.S. (2018) Rapid and Dynamic Alternative Splicing Impacts the Arabidopsis Cold Response Transcriptome. Plant Cell, 30, 1424–1444.
-   Guo,W., Tzioutziou,N., Stephen,G., Milne,I., Calixto,C., Waugh,R., Brown,J.W., and Zhang,R. (2019) 3D RNA-seq - a powerful and flexible tool for rapid and accurate differential expression and alternative splicing analysis of RNA-seq data for biologists. bioRxiv, 656686. doi: https://doi.org/10.1101/656686.
-   Patro,R., Duggal,G., Love,M.I., Irizarry,R.A., and Kingsford,C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods, 14, 417–419.
-   Ritchie,M.E., Phipson,B., Wu,D., Hu,Y., Law,C.W., Shi,W., and Smyth,G.K. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res, 43, e47.
-   Robinson,M.D., McCarthy,D.J., and Smyth,G.K. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139–40.

<a href='#table-of-contents'>Go back to Table of contents</a>

</div>
