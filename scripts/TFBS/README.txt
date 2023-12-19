#Paolo Bagnaresi  
Aim of the script: estimating coherence among proximal Transcription Factor Binding Sites composition versus expression coherence in ortholog gene pairs for barley pantranscriptome members 
script :TFBS_identity_vs_expression_coherence.R 
usage: 
Run the script ensuring that various dependencies as genome files for all genotypes are available and select gene sets of interest 

The script named "TFBS_identity_vs_expression_coherence.R" estimates coherence among proximal Transcription Factor Binding Sites (TFBS) composition versus expression coherence in ortholog gene pairs for barley pan-transcriptome members. 
The script is a wrapper that ultimately calculates percent identities in TFBS in discrete upstream 2Kb regions for gene pairs consisting in Morex gene sets and their close ortholog counterparts and compares such percent identities to expression coherence for the same genes. 
In the first step, given Morex gene sets of interest, close orthologs are identified for each genotype and the presence and distribution patterns of selected 30 TFBS in discrete upstream 2Kb regions are computed and percent identities scored for each gene pair. 
In the second step, for the same close ortholog pairs, the expression coherence (plus or minus 30% coherence in TPM expression values) is computed. Finally, Pearson correlations are calculated among TFBS similarities and expression coherence. 

#Agostino Fricano  
Aim of the script: plotting results of percent identities in TFBS 
script: Plotting_TFBS_identity_coherence.R 
usage:  
Rscript Plotting_TFBS_identity_coherence.R arg1 arg2 arg3 arg3 
arg 1: excel file containing the percentage of identity of cis-element sequences in upstream and downstream regions
arg 2: excel file containing the percentage of TFBS identities against the percentage of coherent expression in each pair comparison 
arg 3: Pearson correlations (in %) between TFBS identity in 2Kb-upstream region and expression coherence 
arg 4: complete path of the output image file with png extension in which plots are saved                                                                               

The script named "Plotting_TFBS_identity_coherence.py" plots the results obtained analyzing the correlation among TFBS Identity coherence in transcript abundance (script TFBS_identity_vs_expression_coherence.R ). This script is written in R and takes four arguments to run: 1) the output of the script "TFBS_identity_vs_expression_coherence.R" reporting the percentage of identity of cis-element sequences in upstream and downstream regions (Excel file), 2) the output of script "TFBS_identity_vs_expression_coherence.R" reporting the percentage of TFBS identities against the percentage of coherent expression in each pair comparison (Excel file), 3) the output of script "TFBS_identity_vs_expression_coherence.R" reporting the Pearson correlations (in %) between TFBS identity in 2Kb-upstream region and expression coherence for 5 sets of genes showing increasing coefficient of variations (Excel file) and 4) the file path in which the image of combined plots is saved in png format. 