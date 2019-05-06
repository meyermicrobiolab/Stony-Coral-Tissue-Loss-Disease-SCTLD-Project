# Stony-coral-tissue-loss-disease-SCTLD-Project

### Requirements
  The packages used and their versions are listed in the sessionInfo.txt file. Some of the key libraries used in analysis are:
  * [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html)
  * [phyloseq](https://joey711.github.io/phyloseq/) 
  * [vegan](https://cran.r-project.org/package=vegan)
  * [CoDaSeq](https://github.com/ggloor/CoDaSeq)
  
### Description 
This project involves the analysis of microbiomes associated with stony coral tissue loss disease from two sites on the Florida Reef Tract. Included in this repository are the original sequencing reads with primers and adaptors removed (also available in NCBI’s Sequence Read Archive under Bioproject PRJNA521988), the R script to recreate the figures in Meyer et al submitted to Frontiers in Microbiology (pre-print available: https://www.biorxiv.org/content/10.1101/626408v1), and the original metadata, ASV, and taxonomy tables from our analysis. The file of session information provides all of the versions of the packages used in R to complete the analysis in Meyer et al. Please note that if the dada2 analysis is repeated, it will result in a slightly different ASV table every time, although the overall interpretations will remain the same. Analysis includes determining the amplicon sequence variants (ASVs) with dada2, assigning taxonomy with the Silva database v.132, and community composition analysis with a combination of phyloseq, vegan, and CoDaSeq tools. 
