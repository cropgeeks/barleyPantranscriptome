setwd("~/PanTrans/PanRTD.clust.c95r0g1s50/R")
options(stringsAsFactors = F)

##################################################
# read in clustering results and gene coordinates
#
##################################################

# cd-hit clustering results in tsv format
clust <- read.table("c95r0g1s50.tab", header = F, sep="\t", col.names = c( "cluster","size", "ID", "ident"))
clust$panID <- do.call(rbind, strsplit(clust$ID, split="::"))[,1]  # get gene ID for later merging

# panRTD gene ID, coordinates and strand.
# first 3 columns in bed format, 4th panRTD gene ID and 5th strand
pan.bed <- read.table("panID.bed", header = F, sep="\t", col.names = c( "chr","start", "end", "panID", "strand"))
clust <- merge(clust,pan.bed, by="panID", all.x=T, sort=F)



#######################################################
# Step 1 filter for clusters with CN>=2

######################################################

# summarize clusters and number of panRTDs within each cluster
clust.nr <- as.data.frame(table(clust$cluster))
colnames(clust.nr)<- c("clust", "CNV")

# CN>=2
clust.nr <- clust.nr[clust.nr$CNV>1,]
CNV2 <- clust[clust$cluster %in% clust.nr$clust, ]




###############################################
# Step 2 
# select clusters on the chromosome 
# calculate range of the gene cluster
# filter for range within 5M
###############################################

cluster <- clust.nr[,1]

df <- data.frame(cluster=cluster, Freq=clust.nr[,2], range=rep(NA, length(cluster)) )


for (i in 1:length(cluster)){
  
  tmp <- CNV2[CNV2$clust == cluster[i],]
  tan <- table(tmp$chr)   # count number of genes on each chr
  if(sum(tan>0)==1){      # only on the same chromosome, calculate range
    tmp <- tmp[tmp$chr %in% names(which(tan>1)),]
    rg <- range(c(tmp$start.bed, tmp$end))
    df$range[i] <- rg[2]-rg[1]
  }
}

df <- df[!is.na(df$range),]   
df <- df[df$range <5000000, ]


CN2.tandem <- CNV2[CNV2$cluster %in% df$cluster, ]


##############################################
## Step 3 
## remove clusters having overlapping genes from same strand 
# which may result in false CNV (+1) if in the same cluster 
############################################
library(GenomicRanges)

plus <- CN2.tandem[CN2.tandem$strand %in% "+" ,]
plus.gr <- GRanges(Rle(plus$chr),IRanges(plus$start+1, plus$end,names=plus$panID ))
plus$group <- subjectHits(findOverlaps(plus.gr , reduce(plus.gr)))
plus.rm <- plus[plus$group %in% which(table(plus$group)>1), ] 


minus <- CN2.tandem[CN2.tandem$strand %in% "-" ,]
minus.gr <- GRanges(Rle(minus$chr),IRanges(minus$start+1, minus$end,names=minus$panID ))
minus$group <- subjectHits(findOverlaps(minus.gr , reduce(minus.gr)))
minus.rm <- minus[minus$group %in% which(table(minus$group)>1), ]  

rem <- unique(c(plus.rm$cluster, minus.rm$cluster))

CN2.tandem <- CN2.tandem[!CN2.tandem$cluster %in% rem , ]




###########################################
# merge cluster info in

###########################################


joint <- merge(CN2.tandem, df, by.x= "cluster", by.y="cluster", all.x=T, sort=F)

write.table(joint, "CN2.tandem.joint.tsv", sep="\t", col.names = T, row.names = F, quote = F)


q("no")























