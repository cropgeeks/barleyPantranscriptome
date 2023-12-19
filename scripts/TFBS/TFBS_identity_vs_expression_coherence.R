
#TFBS_identity_vs_Expression_Coherence
#script for estimation of coherence in proximal CIS-element 2kB region and expression coherence among barley pantranscriptome members

##---required libraries

library("Biostrings")
library("reshape2")
library("seqPattern")

library("ggplot2")
library("writexl")
library("stringr")

###-------------------------
###---custom functions -----
###--------------------------

PF_fas_nodata_no_suffx <- function(DSS, subset_ID , nomefile)
  ### writes fasta file, first arg = data.frame, second arg = ID subset, 
{
   DSS=DSS[names(DSS) %in% subset_ID]
  nomefile = nomefile 
  writeXStringSet(x=DSS, filepath=nomefile,  append=FALSE, format="fasta")
}


PF_MORETT_MOREX_v2_to_PANBART_GENE= function(Id_Morex_V2_QUERY_subset) 
  ### returns IDs of  PANBART gene corresponding to MOREX V2 genes   
{
  return (MORETT_95_95_40[MORETT_95_95_40$Morex_V2_ID   %in%   Id_Morex_V2_QUERY_subset,]$PANBART_GENE  )
}


PF_MORETT_MOREX_v2_to_PANBART_TRANS= function(Id_Morex_V2_QUERY_subset) 
  ###  returns IDs of  PANBART transcripts corresponding to  MOREX V2 genes   
{
  return (MORETT_95_95_40[MORETT_95_95_40$Morex_V2_ID   %in%   Id_Morex_V2_QUERY_subset,]$PANBART_TRANS  )
}

###--------load required data

OK_GsRTD_and_PanBaRT20_match <- read.delim("GsRTD_and_PanBaRT20_match.tsv", header=1)
PanBaRT20 <- read.delim("PanBaRT20.gtf", header=FALSE, comment.char="#")
Horvu_MOREX <- read.delim("Horvu_MOREX.gff", header=FALSE)
AA_MOREX_transfeat_coding_nuc =readAAStringSet('htsm_rtd_Morex_ISO_Morex_RNA_merged_transfeat_coding_pep.fa' )

##loading morex V3 proteins
AA_MOREX_V3_PEP =readAAStringSet('htsm_rtd_Morex_ISO_Morex_RNA_merged_transfeat_coding_pep.fa' )
 
 
##loading morex V2 genome
DSS_MOREX_V2 =readDNAStringSet('Horvu_MOREX_pseudomolecules_v2.fasta.gz' )

##loading morex V2 gff3
MOREX_V2_gff <- read.delim("Horvu_MOREX.gff", header=FALSE)

##--- get gene subset filtering for AHRD quality= 3  and no TE association
ONLYGENES_MOREX_V2_gff= MOREX_V2_gff [MOREX_V2_gff$V3 =="gene",]       
ID_ONLYGENES_MOREX_V2_gff=ONLYGENES_MOREX_V2_gff
ID_ONLYGENES_MOREX_V2_gff$ID=str_replace(ONLYGENES_MOREX_V2_gff$V9,"ID=(.+?);.+$","\\1")
ID_ONLYGENES_MOREX_V2_gff=ID_ONLYGENES_MOREX_V2_gff[, c(10,1:9)]
ID_ONLYGENES_MOREX_V2_gff$AHRD_qual= str_replace(ONLYGENES_MOREX_V2_gff$V9,"^.+?AHRD Quality Score=(.);.+$","\\1")
prelim_te= str_replace(ONLYGENES_MOREX_V2_gff$V9,"^.+?TE-related=(.)?;.+$","\\1")
ID_ONLYGENES_MOREX_V2_gff$TE_related= str_replace(prelim_te,"^.+?TE-related=(.)$","\\1")
ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0 =ID_ONLYGENES_MOREX_V2_gff[ID_ONLYGENES_MOREX_V2_gff$AHRD_qual==3 & ID_ONLYGENES_MOREX_V2_gff$TE_related ==0,]


##---  get CDS NO TE and AHRD qual= 3 
CDS_SOLO_LONGEST_isof_MOREX_V2_HC_CDS =readDNAStringSet('Horvu_MOREX.cds.fasta')  
names(CDS_SOLO_LONGEST_isof_MOREX_V2_HC_CDS)=str_replace(names(CDS_SOLO_LONGEST_isof_MOREX_V2_HC_CDS),"(^.+?)\\..+$","\\1")
CDS_SOLO_LONGEST_isof_MOREX_V2_HC_CDS_AHRD_qual_3__TE_0  =CDS_SOLO_LONGEST_isof_MOREX_V2_HC_CDS[ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0$ID]


##----write to protein  after translation Morex V2 CDS
writeXStringSet (translate(CDS_SOLO_LONGEST_isof_MOREX_V2_HC_CDS_AHRD_qual_3__TE_0,if.fuzzy.codon = "solve" ), "MOREX_V2_AHRD_qual_3__TE_0_PROTEIN.fa")



#####--------------------------------------------------
#####----------get one to one correspondence between  Morex V2 and Morex V3 proteins (data frame named "MORETTA")
#####-----------via best reciprocal hits----------------------------------
 
query_file_name= "/MOREX_V2_AHRD_qual_3__TE_0_PROTEIN.fa"
nome_DATABASE_file ="htsm_rtd_Morex_ISO_Morex_RNA_merged_transfeat_coding_pep.fa" 

MORETTA <- orthologs(query_file     = query_file_name,
                            subject_files     = nome_DATABASE_file,
                            seq_type          = "protein", 
                            ortho_detection   = "RBH",
                            #detailed_output   = TRUE, 
                            eval = "1E-5",
                            clean_folders     = FALSE,
                            comp_cores        = 16,
                            quiet=F)


DF_MORETTA=data.frame(MORETTA)


MORETT_95_95_40 =  DF_MORETTA[DF_MORETTA$evalue <1e-40 &
                             DF_MORETTA$perc_identity >=95 &
                             DF_MORETTA$qcov >= 95  ,]
 
#trim dots  and specify column names
MORETT_95_95_40$ID_trim=str_replace(MORETT_95_95_40$subject_id,"(^.+?)\\..+$","\\1")
names(MORETT_95_95_40 )[1] ="Morex_V2_ID"
names(MORETT_95_95_40 )[2] ="Morex_V3_ID"
names(MORETT_95_95_40 )[22] ="Morex_V3_IDTRIM"

##add PanBartgene and panbarttranscr columns
DF_TEMP_MOR =OK_GsRTD_and_PanBaRT20_match[OK_GsRTD_and_PanBaRT20_match$GsRTD_transcript_id %in%  MORETT_95_95_40$Morex_V3_ID ,] 


MORETT_95_95_40= merge(x =MORETT_95_95_40, y= DF_TEMP_MOR , 
                       by.x="Morex_V3_ID", by.y ="GsRTD_transcript_id",
                       all.x=T)

 
####################################################################################
###---------LOAD TFBS patterns (cis elements consensus) as Biostrings DNAStringSet
#####################################################################################

 RE_PNAS_30 = c("CNGTTR","GKTWGTTR","GKTWGGTR","AGATAT","AGATWCG","GRWAAW","CGCCGCC","RCCGAC","MTCGTA","CAATCA","AATNATT","CGTAC","TTGACY","CAGCT","WGATAR","WATNATW","AAAGY","GGNCCC","TACGTA","GACGTC","CACGTG","TGASTCA","CATGCA","TGTCTC","TGTCGG","CACCTG","TCCGGA","CCTAGG","AACCTA","TAACARA") 
 
 names(RE_PNAS_30 )=   c("MYB-R2R3-MBSI","MYB-R2R3-MBSII","MYB-R2R3-MBSIIG","MYB-related","MYB-GARP-B-ARR","trihelix","AP2-ERF","AP2-DREB","AP2-TOE2","HD-WOX13","HD-ZIP","SPL","WRKY","C2H2-ZF","C2C2-ZF-GATA","C2C2-ZF-YABBY","C2C2-ZF-DOF","TCP","bZIP-A-box","bZIP-C-box","bZIP-G-box","bZIP-GCN4","B3-LEC-AB3","B3-AuxRE","B3-AuxRE2","B3-RAV","LOB","STY1","MRECHS","GARE")
 
 DSS_SEQS_RE_PNAS_30= DNAStringSet( RE_PNAS_30 )
 tag_pattern_set="RE_PNAS_30"
 
 
 
 
 ##############--------------------------------------------------#######
 ##############-----load DATABASE genome-------------------------#######
 ##############--------------------------------------------------#######
 
 genome_DB_files= list.files('genome_V1_folder',pattern = "fasta.gz$", full.names=T)
 genome_DB_gff_files =list.files('genome_V1_gff_folder',pattern="gff$", recursive=T,full.names=T) 
 
tag_QUERY_genotype="MOREX_V2"

#----selection of gene datset 
TRUE_NODES=MORETT_95_95_40$Morex_V2_ID

 DF_ALL_CFR=data.frame()
 str(DF_ALL_CFR)
 
  ###--  ----------------------------------------------
  ### ----------------PATTERN TFBS choice----------------
  ### -----------------------------------------------

 CHOSEN_PATTERNS=DSS_SEQS_RE_PNAS_30

dir_x_ortol_pairs ='directory_for_storage_ortholog_pairs' 

dir.create(dir_x_ortol_pairs, recursive = T)

 ####-------------blast best reciprocal hit     (BRH)                          
  
identity_threshold= 99
query_coverage_HSP_threshold =100
  
  genome_DB_files= list.files('genome_file_folder', pattern="fasta", full.names = T)
  genome_DB_gff_files =list.files('genome_gff_folder', pattern="gff3", full.names = T)

##########################################################################################
##########################################################################################
### part 1   - loop on all genomes for calulation of TFBS percent identity vs Morex  #####
##########################################################################################
##########################################################################################  
 
   for(genome_name in   ALL_GENOMES_DB [-c(16)])   #loop over all genomes apart from Morex
 {
 DSS_DATABASE= NULL
   tag_DB_genotype =genome_name
   cat(" \n loading GENOME : \n  ", basename(genome_DB_files[grep(genome_name, genome_DB_files)]),"\n")
   DSS_DATABASE = readDNAStringSet(genome_DB_files[grep(genome_name, genome_DB_files)] )
   DATABASE_gff= NULL
   cat(" \n loading GFF :  \n ", basename(genome_DB_gff_files[grep(genome_name, genome_DB_gff_files)]),"\n")
   DATABASE_gff <- read.delim(genome_DB_gff_files[ grep(genome_name, genome_DB_gff_files)] , header=FALSE) 
   
   #------create new coluns as ID
   
   ONLYGENES_DATABASE_gff= DATABASE_gff [DATABASE_gff$V3 =="gene",]  
   ID_ONLYGENES_DATABASE_gff=ONLYGENES_DATABASE_gff
   ID_ONLYGENES_DATABASE_gff$ID=str_replace(ONLYGENES_DATABASE_gff$V9,"ID=(.+?);.+$","\\1")
   ID_ONLYGENES_DATABASE_gff=ID_ONLYGENES_DATABASE_gff[, c(10,1:9)]
   ID_ONLYGENES_DATABASE_gff$AHRD_qual= str_replace(ONLYGENES_DATABASE_gff$V9,"^.+?AHRD Quality Score=(.);.+$","\\1")
   prelim_te= str_replace(ONLYGENES_DATABASE_gff$V9,"^.+?TE-related=(.)?;.+$","\\1")
   ID_ONLYGENES_DATABASE_gff$TE_related= str_replace(prelim_te,"^.+?TE-related=(.)$","\\1")
   ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0 =ID_ONLYGENES_DATABASE_gff[ID_ONLYGENES_DATABASE_gff$AHRD_qual==3 & ID_ONLYGENES_DATABASE_gff$TE_related ==0,]
 
  
   CDS_fasta_files = list.files("folder_for_CDS_fasta_files" ,pattern="cds.fasta$", recursive=T,full.names=T)  
   CDS_DATABASE_HC_CDS = readDNAStringSet(CDS_fasta_files[grep(genome_name, CDS_fasta_files)])
  
   #trim .1 
   names(CDS_DATABASE_HC_CDS) =str_replace(names(CDS_DATABASE_HC_CDS),"(^.+?)\\..+$","\\1")

   CDS_DATABASE_HC_CDS_AHRD_qual_3__TE_0  =CDS_DATABASE_HC_CDS[ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0$ID]

 
query_file_name= paste0("Morex_v2_in_MORETT_num_",length(TRUE_NODES),".fa")
PF_fas_nodata_no_suffx(DSS= CDS_SOLO_LONGEST_isof_MOREX_V2_HC_CDS_AHRD_qual_3__TE_0, subset_ID = TRUE_NODES, nomefile = query_file_name)
   
   
   ###
   ####----------------------database file---------------------------------
   ####----------------------database file---------------------------------
   ####----------------------database file---------------------------------
   
    
   nome_DATABASE_file=   paste0("/media/paolo/SEAGATE/paolo/AAA_NUOVI_2K20/AAA_PANGEN_BARLEY/AA_NUOVI_DATI_LUG_23/FOLDER_CORSE/FASTA_files/DATABASE_FASTA_CDS_SOLO_LONGEST_isof_HC_CDS_AHRD_qual_3__TE_0.fa" )
   writeXStringSet(CDS_DATABASE_HC_CDS_AHRD_qual_3__TE_0,filepath = nome_DATABASE_file)
   
   
   
   ###   use old ortholog pairs ? 
   if (use_dir_old_orto_pairs==F)
   {  
    ####----------------------database file---------------------------------

     
     library(orthologr)
     
     LS_DF_MERGE_RBH=list()
     
     RBH_pao <- orthologs(query_file     = query_file_name,
                          subject_files     = nome_DATABASE_file,
                          seq_type          = "cds", 
                          ortho_detection   = "RBH",
                          #detailed_output   = TRUE, 
                          eval = "1E-2",
                          clean_folders     = FALSE,
                          comp_cores        = 13,
                          quiet=F)
     
     
     DF_RBH_pao=data.frame(RBH_pao)
     
##------------hortolog pairs selection
     
DF_RBH_pao=data.frame(RBH_pao)
DF_RBH_pao = DF_RBH_pao[grep("_Un",DF_RBH_pao$query_id, invert = T),]   #elimino unidentified in query
DF_RBH_pao_OK=  DF_RBH_pao[grep("_Un",DF_RBH_pao$subject_id, invert = T),] #elimino unidentified in subj
DF_RBH_pao_OK_FINALE_BIS = DF_RBH_pao_OK[DF_RBH_pao_OK$perc_identity >=  identity_threshold & 
                                           DF_RBH_pao_OK$qcovhsp >= query_coverage_HSP_threshold  & 
                                           (DF_RBH_pao_OK$q_end - DF_RBH_pao_OK$q_start +1) ==DF_RBH_pao_OK$q_len &
                                           (DF_RBH_pao_OK$s_end - DF_RBH_pao_OK$s_start +1) == DF_RBH_pao_OK$s_len ,]
     

##------RBH storage in list-----
 
 LS_DF_RBH_pao_OK_FINALE_BIS[[paste0("TUTTI_MOREX_V2_IN_MORETTA_VS_",tag_DB_genotype)]]= DF_RBH_pao_OK_FINALE_BIS
     save(LS_DF_RBH_pao_OK_FINALE_BIS, file= paste0('/media/paolo/SEAGATE/paolo/AAA_NUOVI_2K20/AAA_PANGEN_BARLEY/AA_NUOVI_DATI_LUG_23/FOLDER_CORSE/RBH_OGG_R/',gene_set,'.RData' ))
     
   }   
   

   ##-----------------------------------------------------------------------------------------------------------------------
   ##-----------------------ordering seqs  ------------
   ##-----------------------------------------------------------------------------------------------------------------------
   
   vec_order_seqs=NULL
   
   for (hh in 1:nrow(  DF_RBH_pao_OK_FINALE_BIS))
     # hh=1
   {vec_order_seqs  =c(vec_order_seqs, paste0( DF_RBH_pao_OK_FINALE_BIS[hh,1:2])) }
   
   TRUE_NODES_QUERY =    DF_RBH_pao_OK_FINALE_BIS$query_id
   TRUE_NODES_DATABASE =   DF_RBH_pao_OK_FINALE_BIS$subject_id
   
   ####-----------------------------------------------------------------
   ####----cration of DNAStringSets for retrieval of upstream regions                
   ####-----------------------------------------------------------------
   
   NT_UPSTREAM=2000
   NT_DOWNSTREAM=1000
 
   ##------------------------------------------------------------------------------------  
   ##  -------MOREX FORWARD--------------------------------------------------------------------    
   ##  -------MOREX FORWARD--------------------------------------------------------------------  
   ##------------------------------------------------------------------------------------
   
SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW = ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0[ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0$ID %in%  DF_RBH_pao_OK_FINALE_BIS$query_id & ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0$V7=="+",]
  
    
   system.time (
     {   
      DF_sint_seqs_prom_FORW=NULL
    l_prealloc =length(SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW$ID)
       
  DF_sint_seqs_prom_FORW=data.frame(ID_gene =  character(l_prealloc), chrom= character( l_prealloc), gene_start=integer ( l_prealloc),
                                         gene_end= integer  ( l_prealloc), seq= character( l_prealloc),strand=character ( l_prealloc))
       
       
       for (step in 1: l_prealloc)
       { 
        DF_sint_seqs_prom_FORW[step,1 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW$ID[step] 
        DF_sint_seqs_prom_FORW[step,2 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW$V1[step]
        DF_sint_seqs_prom_FORW[step,3 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW$V4 [step] 
        DF_sint_seqs_prom_FORW[step,4 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW$V5[step] 
      DF_sint_seqs_prom_FORW[step,5 ]= as.character( DNAStringSet(DSS_MOREX_V2[ SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW$V1[step]],
                                                                     start = SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW[step,]$V4 -NT_UPSTREAM, 
                                                                     width=  NT_UPSTREAM  + NT_DOWNSTREAM))
         DF_sint_seqs_prom_FORW[step,6]  =SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_FORW$V7[step] 
       }}
     
   )
   
   
   
   
   DSS_seqs_FORW =DNAStringSet(DF_sint_seqs_prom_FORW$seq)
   
   names(DSS_seqs_FORW)=paste0( DF_sint_seqs_prom_FORW$ID_gene,"__gene_START_",DF_sint_seqs_prom_FORW$gene_start, "__UP_",NT_UPSTREAM ,"__DOWN_",NT_DOWNSTREAM )
   DSS_seqs_FORW
   
   
   ##--MOREX_REVCOMPL-------------------------------------
   ##--MOREX_REVCOMPL-------------------------------------
   
  SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL= ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0[ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0$ID %in% DF_RBH_pao_OK_FINALE_BIS$query_id & ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0$V7=="-",]
   
   # str(SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL)  #6046 more v2  4ago23
   
   
   system.time({  
     
     DF_sint_seqs_prom_REVCOMPL=NULL
     ##inizializzo DF revcompl
     l_prealloc =length(SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL$ID)
     
     DF_sint_seqs_prom_REVCOMPL=data.frame(ID_gene =  character(l_prealloc),chrom=character( l_prealloc),gene_start=integer ( l_prealloc),
                                           gene_end=      integer  ( l_prealloc),  seq= character( l_prealloc),strand=    character ( l_prealloc))
     
     
     for (step in 1: l_prealloc)
     { 
       DF_sint_seqs_prom_REVCOMPL[step,1 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL$ID[step] 
       DF_sint_seqs_prom_REVCOMPL[step,2 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL$V1[step]
       DF_sint_seqs_prom_REVCOMPL[step,3 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL$V4 [step] 
       DF_sint_seqs_prom_REVCOMPL[step,4 ]= SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL$V5[step] 
       DF_sint_seqs_prom_REVCOMPL[step,5 ]= as.character( reverseComplement(DNAStringSet(DSS_MOREX_V2[SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL$V1[step]],
                                                                                         start = SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL[step,]$V5 -NT_DOWNSTREAM, 
                                                                                         width=  NT_UPSTREAM  + NT_DOWNSTREAM)))
       
       DF_sint_seqs_prom_REVCOMPL[step,6]  =SELEZ_ID_ONLYGENES_MOREX_V2_gff_AHRD_qual_3__TE_0_REVCOMPL$V7[step] 
     }}
     
   )
   
 

   DSS_seqs_REVCOMPL=DNAStringSet(   DF_sint_seqs_prom_REVCOMPL$seq)
   
   names(DSS_seqs_REVCOMPL)=paste0( DF_sint_seqs_prom_REVCOMPL$ID_gene,"__gene_START_",DF_sint_seqs_prom_REVCOMPL$gene_start, "__UP_",NT_UPSTREAM ,"__DOWN_",NT_DOWNSTREAM )
   DSS_seqs_REVCOMPL
   
 ####-------------------  merging forward and reverse sequences              
   
DSS_MERGE_FORW_REV = append(DSS_seqs_FORW,DSS_seqs_REVCOMPL)
   
####-------------------  merging forward and reverse sequences   


################################################
###############  DATABASE genotypes   ##########
################################################

#  -----DATABASE seqs FORWARD------------------------------------------------------------------------------                  
#  -----DATABASE seqs FORWARD------------------------------------------------------------------------------                  

   
SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW = ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0[ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0$ID %in% DF_RBH_pao_OK_FINALE_BIS$subject_id & ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0$V7=="+",]
   str( SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW) #aksh 6316 4ago23
   
   system.time (
     {   
       DF_sint_seqs_DATABASE_prom_FORW=NULL
       l_prealloc =length( SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW$ID)
       
       DF_sint_seqs_DATABASE_prom_FORW =data.frame(ID_gene =  character(l_prealloc), chrom= character( l_prealloc), gene_start=integer ( l_prealloc),
                                                   gene_end= integer  ( l_prealloc), seq= character( l_prealloc),strand=character ( l_prealloc))
       
       
       for (step in 1: l_prealloc)
       { 
         DF_sint_seqs_DATABASE_prom_FORW[step,1 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW$ID[step] 
         DF_sint_seqs_DATABASE_prom_FORW[step,2 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW$V1[step]
         DF_sint_seqs_DATABASE_prom_FORW[step,3 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW$V4 [step] 
         DF_sint_seqs_DATABASE_prom_FORW[step,4 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW$V5 [step] 
         DF_sint_seqs_DATABASE_prom_FORW[step,5 ]= as.character( DNAStringSet(DSS_DATABASE[ SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW$V1[step]],
                                                                              start = SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW[step,]$V4 -NT_UPSTREAM, 
                                                                              width=  NT_UPSTREAM  + NT_DOWNSTREAM))
         DF_sint_seqs_DATABASE_prom_FORW[step,6]  =SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_FORW$V7[step] 
       }}
     )
   
  
   DSS_seqs_DATABASE_FORW  =DNAStringSet(DF_sint_seqs_DATABASE_prom_FORW$seq)
   
names(DSS_seqs_DATABASE_FORW)=paste0( DF_sint_seqs_DATABASE_prom_FORW$ID_gene,"__gene_START_",DF_sint_seqs_DATABASE_prom_FORW$gene_start, "__UP_",NT_UPSTREAM ,"__DOWN_",NT_DOWNSTREAM )
   
   

##----------------------------------------------------------
##---genotype DATABASE REVCOMPL-------------------------------------
##---------------------------------------------------------

SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL = ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0[ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0$ID %in% DF_RBH_pao_OK_FINALE_BIS$subject_id & ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0$V7=="-",]
   str(  SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL ) 
   
   system.time({  
     
     DF_sint_seqs_DATABASE_prom_REVCOMPL=NULL
     ##inizializzo DF revcompl
     l_prealloc =length(SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL$ID)
     
     DF_sint_seqs_DATABASE_prom_REVCOMPL=data.frame(ID_gene =  character(l_prealloc),chrom=character( l_prealloc),gene_start=integer ( l_prealloc),
                                                    gene_end=      integer  ( l_prealloc),  seq= character( l_prealloc),strand=    character ( l_prealloc))
     
     
     for (step in 1: l_prealloc)
     { 
       DF_sint_seqs_DATABASE_prom_REVCOMPL[step,1 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL$ID[step] 
       DF_sint_seqs_DATABASE_prom_REVCOMPL[step,2 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL$V1[step]
       DF_sint_seqs_DATABASE_prom_REVCOMPL[step,3 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL$V4[step] 
       DF_sint_seqs_DATABASE_prom_REVCOMPL[step,4 ]= SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL$V5[step] 
       DF_sint_seqs_DATABASE_prom_REVCOMPL[step,5 ]= as.character( reverseComplement(DNAStringSet(DSS_DATABASE[ SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL$V1[step]],
                                                                                                  start = SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL[step,]$V5 -NT_DOWNSTREAM, 
                                                                                                  width=  NT_UPSTREAM  + NT_DOWNSTREAM)))
       
       DF_sint_seqs_DATABASE_prom_REVCOMPL[step,6]  =SELEZ_ID_ONLYGENES_DATABASE_gff_AHRD_qual_3__TE_0_REVCOMPL$V7[step] 
     }}
     
   )
   
   str(DF_sint_seqs_DATABASE_prom_REVCOMPL)
   DSS_seqs_DATABASE_REVCOMPL=DNAStringSet(   DF_sint_seqs_DATABASE_prom_REVCOMPL$seq)
   
   names( DSS_seqs_DATABASE_REVCOMPL)=paste0( DF_sint_seqs_DATABASE_prom_REVCOMPL$ID_gene,"__gene_START_",DF_sint_seqs_DATABASE_prom_REVCOMPL$gene_start, "__UP_",NT_UPSTREAM ,"__DOWN_",NT_DOWNSTREAM )
   DSS_seqs_DATABASE_REVCOMPL
   
   
   
   
   
    ##----MERGING DATABASE FORW and REVCOMPL                                                                                                             
   
   DSS_MERGE_DATABASE_FORW_REV = append(DSS_seqs_DATABASE_FORW, DSS_seqs_DATABASE_REVCOMPL)
   DSS_MERGE_DATABASE_FORW_REV
   

   ###----merging Morex and Database seqs
 
   
   DSS_MERGEALL=NULL

   DSS_MERGEALL =append( DSS_MERGE_FORW_REV, DSS_MERGE_DATABASE_FORW_REV)
   DSS_MERGEALL_idtrim=DSS_MERGEALL
   names(DSS_MERGEALL_idtrim) =str_replace(names(DSS_MERGEALL),"(^.+?)__.+$","\\1")
  
##---sequence  pairing loop
   
   DSS_testorder=DNAStringSet()
   for (a in vec_order_seqs)
   {  #print(a)
     DSS_testorder = append( DSS_testorder,DSS_MERGEALL_idtrim[a]) }
   
   #some cleaning 
   rm(DSS_MERGEALL_idtrim)
   rm(DSS_MERGEALL)
   gc()

  ###----- ----filter for seqs containg "Ns"--------
  CON_N =grep("N", DSS_testorder)
   
   if (!identical(  CON_N , integer(0)))
   {  
      to_be_cancel=NULL
   
   for (num in CON_N)
   {
     to_be_cancel= c(to_be_cancel, num)
     if((num %% 2) == 0) { to_be_cancel= c(to_be_cancel, num -1)}
     else { to_be_cancel= c(to_be_cancel, num + 1)}
   }
   s_to_be_cancel=to_be_cancel[order(to_be_cancel, decreasing =F)]
   
   #####-----------------------------------------------------------
   
   
DSS_testorder_NO_N= DSS_testorder[-c(s_to_be_cancel)]
DSS_testorder_orig_pre_NON =DSS_testorder  # keep track of all sequences
 DSS_testorder=DSS_testorder_NO_N
     }


   DSS_testorder_NO_N = DSS_testorder 
   ###-----------------------------------------------
   ###--write down paired sequences after N filtering
###--------------------------------------------------
   
   DF_RBH_per_scrittura_no_Ns =DF_RBH_pao_OK_FINALE_BIS[DF_RBH_pao_OK_FINALE_BIS$query_id %in% names(DSS_testorder_NO_N) &
                                                          DF_RBH_pao_OK_FINALE_BIS$subject_id %in% names(DSS_testorder_NO_N) , ]
   
   
   write_xlsx(DF_RBH_per_scrittura_no_Ns, path = "ARBITRARY_PATH")
   
   ### prepare for TFBS identification 
  
   library(reshape2)
   library(seqPattern)
   
   tag_QUERY_genotype="MOREX_V2"
   DSS_seqs= DSS_testorder_NO_N
    write_or_not_xls="yes"

  
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   LS_DF_GLOBAL= list()
   HYPER_MELT_DF_STATS_CONSIST_GLOBAL= data.frame() 
   CAST_DF_STATS_CONSIST_XLS_GLOBAL=data.frame()
   
  bb=NULL
   for (bb in LETTERS[1:6] )
     
   {
     if(bb=="A") { DSS_seqs_STATS= DNAStringSet(DSS_seqs,start=1,width=500)  
     regione = paste0("A_", length(DSS_seqs)/2,"__", tag_DB_genotype) ## ###  up streamcoding 1000bp 
     regione_short="A"
     }
     
     if(bb=="B")  {DSS_seqs_STATS= DNAStringSet(DSS_seqs,start=501,width=500)  
     regione= paste0("B_",length(DSS_seqs)/2, "__", tag_DB_genotype) ## ###  up streamcoding 1000bp 
     regione_short="B"
     }
     
     if(bb=="C")  {DSS_seqs_STATS= DNAStringSet(DSS_seqs,start=1001,width=500)  
     regione= paste0("C_",length(DSS_seqs)/2, "__", tag_DB_genotype) ## ###  up streamcoding 1000bp   
     regione_short="C"
     }
     
     if(bb=="D")   {DSS_seqs_STATS= DNAStringSet(DSS_seqs,start=1501,width=500)
     regione= paste0("D_",length(DSS_seqs)/2, "__", tag_DB_genotype) 
     regione_short="D"
     }
     if(bb=="E")   {DSS_seqs_STATS= DNAStringSet(DSS_seqs,start=2001,width=500)
     regione= paste0("E_",length(DSS_seqs)/2, "__", tag_DB_genotype) 
     regione_short="E"
     }
     if(bb=="F")   {DSS_seqs_STATS= DNAStringSet(DSS_seqs,start=2501,width=500)
     regione= paste0("F_",length(DSS_seqs)/2, "__", tag_DB_genotype) 
     regione_short="F"
     }
     
     cat("\n\n doing ", regione, "\n\n")
    
     DF_STATS_CONSIST_GLOBAL=NULL
     DF_STATS_CONSIST_XLS_GLOBAL=NULL
     
     DF_stats_CONSIST=NULL
     DF_stats_CONSIST_XLS=NULL
     OCC_list_QUERY=NULL
     OCC_list_DATABASE = NULL

 
###---identify TFBS occurrences and positioning by seqPattern function  "getPatternOccurrenceList"

OCC_list_QUERY = getPatternOccurrenceList(DSS_seqs_STATS[seq(from = 1, to =length(DSS_seqs_STATS),by=2)] ,patterns=  CHOSEN_PATTERNS   ,  useMulticore = T, nrCores = num_cores_for_seqpattern) #cores erano 15 !!  10mar21
     OCC_list_QUERY_NONULL =OCC_list_QUERY [!sapply(OCC_list_QUERY , is.null)]
### database genotype     
OCC_list_DATABASE = getPatternOccurrenceList(DSS_seqs_STATS[seq(from = 0, to =length(DSS_seqs_STATS),by=2)] ,patterns=  CHOSEN_PATTERNS   ,  useMulticore = T, nrCores = num_cores_for_seqpattern)
     OCC_list_DATABASE_NONULL =OCC_list_DATABASE [!sapply(OCC_list_DATABASE , is.null)]
     
 ### ----loop over the 30 TFBS
     for (patt in   as.character(CHOSEN_PATTERNS ) )
     {
       pasted_hits_q = paste0(OCC_list_QUERY_NONULL[[patt]][,1],"-",OCC_list_QUERY_NONULL[[patt]][,2])
       pasted_hits_d = paste0(OCC_list_DATABASE_NONULL[[patt]][,1],"-",OCC_list_DATABASE_NONULL[[patt]][,2])

       pasted_hits_q_num = as.integer(paste0(OCC_list_QUERY_NONULL[[patt]][,1],OCC_list_QUERY_NONULL[[patt]][,2]))
       pasted_hits_d_num = as.integer(paste0(OCC_list_DATABASE_NONULL[[patt]][,1],OCC_list_DATABASE_NONULL[[patt]][,2]))
       

      
 ###--datafrane for summarizing TFBS hits      

       DF_stats_CONSIST  = data.frame(    
         TOTAL_hits_in_QUERY= length( OCC_list_QUERY_NONULL[[patt]][,1]) ,
         TOTAL_hits_in_DB= length( OCC_list_DATABASE_NONULL[[patt]][,1]),
         REGION=  regione_short,
         CIS_ELEM_name =  paste0(names(CHOSEN_PATTERNS[CHOSEN_PATTERNS==patt ])," - ",as.character(patt)), 
    
         N_SET_perf_hits= length(intersect(pasted_hits_q , pasted_hits_d)),
         hits_missing_in_query= length(setdiff(pasted_hits_q , pasted_hits_d)),
         hits_missing_in_db = length(setdiff(pasted_hits_d , pasted_hits_q))
       )
      
       DF_STATS_CONSIST_GLOBAL= rbind(DF_STATS_CONSIST_GLOBAL, DF_stats_CONSIST)
       
     }
     
  
     
     
     MELT_DF_STATS_CONSIST_GLOBAL <- melt(DF_STATS_CONSIST_GLOBAL,
                                          id.vars = c("REGION","CIS_ELEM_name"),
                                          variable.name = "Type",
                                          value.name = "Value")
     
     LS_DF_GLOBAL[[regione]] = MELT_DF_STATS_CONSIST_GLOBAL
     

     #....................................#PLOT CIS STATS.............................................
     #....................................#PLOT CIS STATS.............................................
     #....................................#PLOT CIS STATS.............................................
     
     LS_DF_GLOBAL[[regione]]
     HYPER_MELT_DF_STATS_CONSIST_GLOBAL <- do.call("rbind", LS_DF_GLOBAL)
   
       gc()
     
 

     
   }

 
   ##--- preapare dat frame to get a spreadsheet
   
   uuT=dcast(HYPER_MELT_DF_STATS_CONSIST_GLOBAL,...~Type)

   
   UU_perc= uuT
   
   
   UU_perc$GROUP=nome_GROUP
   UU_perc$CONDITION=gene_set
   UU_perc$GENOTYPE =str_replace(tag_DB_genotype,"^.+_(.+$)","\\1")
   UU_perc$SET_perc_perfect_hits= uuT$N_SET_perf_hits/uuT$TOTAL_hits_in_QUERY *100
   UU_perc$perc_hits_in_query_missing_in_DB = uuT$hits_missing_in_query/uuT$TOTAL_hits_in_QUERY*100
   UU_perc$perc_hits_in_DB_missing_in_query = uuT$hits_missing_in_db/uuT$TOTAL_hits_in_QUERY*100
   UU_perc$PAIRS =length(DSS_seqs)/2
   UU_perc$NHITS_100seqsPAIRS = uuT$N_SET_perf_hits/(length(DSS_seqs)/2)*100
   UU_perc$log = log10(uuT$N_SET_perf_hits/(length(DSS_seqs)/2)*100)
   UU_perc$ln = log(uuT$N_SET_perf_hits/(length(DSS_seqs)/2)*100)
   

   
UU_perc_FINAL= UU_perc[c("REGION","CONDITION","GROUP","GENOTYPE","CIS_ELEM_name", "TOTAL_hits_in_QUERY","TOTAL_hits_in_DB","N_SET_perf_hits", "hits_missing_in_query" ,"hits_missing_in_db","SET_perc_perfect_hits","perc_hits_in_query_missing_in_DB","perc_hits_in_DB_missing_in_query","PAIRS","NHITS_100seqsPAIRS","log","ln" )]
   
   
   #.................................................................................
   ##----------------------------------write xlsx file-----------------
   #...................................................................................
   
file_name_UU_perc_FINAL= paste0(tag_QUERY_genotype,"_VS_",tag_DB_genotype,"__CERC_",gene_set,"_patterns_",tag_pattern_set)
   
folder_all_comparisons= "final_directory_for_writing_%_TFBS_identities"
   

write_xlsx(UU_perc_FINAL, path = paste0( folder_all_comparisons ,file_name_UU_perc_FINAL,  "_seq_pairs_",length(DSS_seqs)/2,"_",PF_ora(),".xlsx"))
 

   
   cat("\n\n\ finished  ---> ", tag_DB_genotype, "\n")
   
   
  ###-------------------------------------------------------------------
  ###---------------write down all analyses ---------------------------
  ###-----------------------------------------------------------------
   
   DF_ALL_CFR= rbind(  DF_ALL_CFR, UU_perc_FINAL)
                  
   
   write_xlsx(DF_ALL_CFR, path = paste0(folder_all_comparisons ,"DATAFRAME_ALL.xlsx"))
 
   }
  
  
##########################################
###     end part 1               #########
########################################## 
  
  
  
#############################################################################
############################################################################# 
###  part 2    expression coherence   calulations                     #######
#############################################################################  
#############################################################################  
  
  
  
AA_MAIN_DF_ALL_CFR   <- read_excel(  DATAFRAME_ALL) #read all data of TFBS pecent idntity, condition
DIR_ORTHOL_PAIRS= # directory where all ortolog pairs are stored
DIR_o_CFR_MOREX_vs= # set  directory for new data
 
  choice_specific_cis_elem= F # consider specifically one or more  cis-elem (TFBS) or all 30
  vector_upstream_region= c("A|B|C|D" )  #selection of one or more upstream regions
 
  
for (upstream_region in vector_upstream_region) 
  { 
    for (percen_coer in c( 30))
    {
      threshold_up  =1+ percen_coer*1/100 
      threshold_down =1-percen_coer*1/100  
      
      cat("\nthreshold_up= ",threshold_up)
      cat("\nthreshold_down= ",threshold_down)
      
      accu_tpm_coeren=NULL
      accu_num_ortholog_pairs=NULL
      accu_percen_coer=NULL
      accu_cfr=NULL
      accu_percen_30cis_REG_CHOSEN=NULL
      DF_temp=NULL
 
      
      for (genome_name_for_All in nomi_genom_adatti_TPM_MORETT)   
      { 
        cat("\n\n  mome_genoma  ",genome_name_for_All ,"\n\n")
 
        if ( genome_name_for_All != "Morex") # avoid Morex 
        {
          #load gene IDs
          Id_Morex_V2_QUERY_subset=NULL
          Id_Morex_V2_QUERY_subset <-read_excel(list.files(DIR_ORTHOL_PAIRS, pattern=paste0( searched_condition ,"_.+?",genome_name_for_All,"_BEST_REC_HITS_" ),full.names=T)[1])$query_id
          length( Id_Morex_V2_QUERY_subset)       
          
          #number ortholog pairs
          num_orthol_pairs = length( Id_Morex_V2_QUERY_subset)
 
          PANBART_GENES_subset= PF_MORETT_MOREX_v2_to_PANBART_GENE( Id_Morex_V2_QUERY_subset) ##get panbart gene IDS
 
          ###.......................................MOREX rowmeans ...................
          ###.......................................MOREX rowmeans ...................
          ###.......................................MOREX rowmeans ...................
          
          MOREX_rowmeans=NULL
          #get Morex rowmeans tpm for gene set of interesin panbart 
          MOREX_rowmeans = rowMeans(  TPM_genes_inMORETT[ PANBART_GENES_subset, grep("Morex",   colnames(TPM_genes_inMORETT))]) 
          
          ### .................COUNTERPART....................
          ###..................COUNTERPART....................
          ###..................COUNTERPART....................
          
          #load counterpart  data to compare to  morex e of same subset
          GENERIC_rowmeans =NULL     
          GENERIC_rowmeans= rowMeans(TPM_genes_inMORETT[PANBART_GENES_subset,
                                                        grep(genome_name_for_All, colnames(TPM_genes_inMORETT))])
   
        ##--------------------------- 
        ##estimate TPM coherence
        ##--------------------------
          
        tpm_coeren=NULL
        tpm_coeren=length(which( MOREX_rowmeans/GENERIC_rowmeans <= threshold_up & MOREX_rowmeans/GENERIC_rowmeans>= threshold_down))
        cat( "\nMOREX_V2 fold expression  vs ", genome_name_for_All, " orthol pairs ", num_orthol_pairs, "within fold expression range ",   threshold_up , " and ",threshold_down, " = " , tpm_coeren, "    -----> percent= ", tpm_coeren/length(MOREX_rowmeans)*100,"\n\n")

          
          if (choice_specific_cis_elem == F)   
          { 
            
            DF_temp =NULL
            DF_temp =  AA_MAIN_DF_ALL_CFR[grepl(upstream_region, AA_MAIN_DF_ALL_CFR$REGION) & 
                                              AA_MAIN_DF_ALL_CFR$CONDITION ==searched_condition   &
                                              AA_MAIN_DF_ALL_CFR$GENOTYPE == genome_name_for_All,]
          }
          
          
          if (choice_specific_cis_elem== T)  
          {
           cis_elem_chosen = "MYB-GARP-B-ARR - AGATWCG" 
            DF_temp =NULL
            DF_temp =  AA_MAIN_DF_ALL_CFR[grepl(upstream_region, AA_MAIN_DF_ALL_CFR$REGION)   & 
                                              AA_MAIN_DF_ALL_CFR$CIS_ELEM_name  == cis_elem_chosen &
                                              AA_MAIN_DF_ALL_CFR$CONDITION ==searched_condition    &
                                              AA_MAIN_DF_ALL_CFR$GENOTYPE == genome_name_for_All,]
            
          }
          
          
          
          #........................create final dataframe  -------
          
          accu_tpm_coeren = c(accu_tpm_coeren, tpm_coeren)
          accu_num_ortholog_pairs=c(accu_num_ortholog_pairs,num_orthol_pairs)
          accu_percen_coer=c(accu_percen_coer,tpm_coeren/length(MOREX_rowmeans)*100)
          accu_cfr=c(accu_cfr,  paste0("Morex__vs__",genome_name_for_All ,"___",num_orthol_pairs,"_pairs"))
          accu_percen_30cis_REG_CHOSEN=c(accu_percen_30cis_REG_CHOSEN, mean(DF_temp$percent_INTERSECT_perf_hits))  
          
          
        }    
  
      }  
      
     
      
      
      ####   DATAFRAME----------------------------......
      
      cfr_morex_vs_others= data.frame(tpm_coeren=accu_tpm_coeren,
                                     perc_coer=accu_percen_coer,
                                  
                                     cfr= accu_cfr,
                                     orth_pairs= accu_num_ortholog_pairs,
                                     perc_ident_30cis_REGION_CHOSEN= accu_percen_30cis_REG_CHOSEN)
      
      ###---order-plots----------------------------------------
 
      o_cfr_morex_vs_others =    cfr_morex_vs_others[order(    cfr_morex_vs_others$tpm_coeren, decreasing=F),  ]
      #LIST_o_cfr=list()  store data in list
      LIST_o_cfr[[paste0(searched_condition)]]=    o_cfr_morex_vs_others

      write_xlsx( LIST_o_cfr[[searched_condition]] ,paste0( DIR_o_CFR_MOREX_vs, searched_condition ,".xlsx"))
      
  
  
      ######################################################## 
      ####   get   pearson correlations and save ############
      ########################################################
      
      
  cor_classic=round(cor(x = o_cfr_morex_vs_others$perc_coer,  y =   o_cfr_morex_vs_others$perc_ident_30cis_REGION_CHOSEN),3)
      
      ##save cor_classic as necessary
      #LIST_o_cfr=list()
      LIST_o_cfr[[paste0(searched_condition)]]=       o_cfr_morex_vs_others
      names(LIST_o_cfr)

      write_xlsx( LIST_o_cfr[[searched_condition]] ,paste0( DIR_o_CFR_MOREX_vs, searched_condition ,".xlsx"))
      

  
    }
}
  
  
  
 