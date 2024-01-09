# Barley pan-transcriptome paper -- analysis code
 
This repository contains the analysis scripts used in the barley pan-transcriptome paper.

These are organised as follows:
- scripts/geneClustering/Gene.cluster.filter.R: this script is used to filter gene duplication clusters. It reads gene clustering results and gene info (coordinates and strand) as input files. It outputs a table of gene clusters in tandem, cluster range and number of genes in the cluster.
- scripts/MorexAtlas: contains all scripts used for the construction of the Morex gene expression atlas and for hosting it on the EORNA web infrastructure.
- scripts/TFBS: contains the scripts used for analysing the transcription factor binding sites and plotting the results from this.

  There is a separate [README](https://github.com/cropgeeks/barleyPantranscriptome/blob/main/PanBaRT20.md) that details the construction of the PanBaRT20 RTD.
  
  Scripts used in the network construction and analysis are found [here](https://github.com/vanda-marosi/PanBarleyNetworks/tree/5fce40c52965ab6bd4d5c45e1438a5fa1ab36a4f).
  
  Scripts used for the phenology gene analysis are hosted in a separate github repository [here](https://github.com/WCGA-Murdoch/Barley-phenology-2023).
