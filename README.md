# Stony-Coral-Tissue-Loss-Disease-SCTLD-Project

### Requirements
  The packages used and their versions are listed in the sessionInfo.txt file. Some of the key libraries used in analysis are:
  * [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html)
  * [phyloseq](https://joey711.github.io/phyloseq/) 
  * [vegan](https://cran.r-project.org/package=vegan)
  * [CoDaSeq](https://github.com/ggloor/CoDaSeq)
  
### Description 
This is the set of R scripts used for the analysis of microbiomes associated with stony coral tissue loss disease from two sites on the Florida Reef Tract. These scripts will allow the user to recreate the figures in Meyer et al submitted to Frontiers in Microbiology. Sequencing reads with adapters and primers removed are available in this repository as well as in NCBI’s Sequence Read Archive under Bioproject PRJNA521988. Analysis includes determining the amplicon sequence variants (ASVs) with dada2, assigning taxonomy with the Silva database v.132, and community composition analysis with a combination of phyloseq, vegan, and CoDaSeq tools. 
