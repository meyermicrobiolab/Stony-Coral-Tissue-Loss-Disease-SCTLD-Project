library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(vegan)
library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(cowplot)
library(randomcoloR)
library(DESeq2)
library(dplyr)
library(reshape2)
library(tibble)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

###### Quality-filter reads and create Amplicon Sequence Variant tables
# put parsed, adaptors & primers removed, unjoined (R1 and R2 separate) fastq files
# into directory for DADA2 & make sure the full path is updated in the next line:

###adjust path name for testing
path <- "~/Documents/DiseaseOutbreak/DiseaseMicrobiomeAnalysis/FIELDSITES/cutadapt_fieldsites"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_cut.fastq and SAMPLENAME_R2_cut.fastq
# Samplename is everything before the first underscore
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

# Perform filtering and trimming
# Assign the filenames for the filtered fastq.gz files.
# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the Error Rates, it TAKES TIME! do first forward and then reverse
# Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

# visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspecting the dada-class object returned by dada:
dadaFs[[1]]

# Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers) ## The sequences being tabled vary in length.
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
# As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "dada_read_stats.txt",sep="\t",col.names=NA)

# SAVE THIS FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS, adjust name
saveRDS(seqtab.nochim, file="~/Documents/DiseaseOutbreak/DiseaseMicrobiomeAnalysis/FIELDSITES/seqtab.nochim.rds")
# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
# seqtab.nochim <- readRDS("~/Documents/DiseaseOutbreak/DiseaseMicrobiomeAnalysis/FIELDSITES/seqtab.nochim.rds")

# Assign taxonomy
# Make sure the appropriate database is available in the DADA2 directory
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/DiseaseOutbreak/DiseaseMicrobiomeAnalysis/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# FIX the NAs in the taxa table
taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]
write.table(taxon,"Fieldsites_silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "Fieldsites_silva_otu_table.txt",sep="\t",col.names=NA)

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata.
otu <- read.table("Fieldsites_silva_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Fieldsites_silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Fieldsites_metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
#12300 taxa and 65 samples

# remove chloroplasts and mitochondria and Eukaryota
get_taxa_unique(ps, "Family") #597
get_taxa_unique(ps, "Order") #339
get_taxa_unique(ps, "Kingdom") #4
ps2 <- subset_taxa(ps, Family !="Mitochondria")
ps2 <- subset_taxa(ps2, Order !="Chloroplast")
ps2 <- subset_taxa(ps2, Kingdom !="Eukaryota")
ps2 <- subset_taxa(ps2, Kingdom !="NA")
get_taxa_unique(ps2, "Family") #593
get_taxa_unique(ps2, "Order") #336
get_taxa_unique(ps2, "Kingdom") #2
ps2
#11332 taxa and 65 samples

# filtered taxa with phyloseq, now export cleaned otu and taxa tables from phyloseq
otu = as(otu_table(ps2), "matrix")
taxon = as(tax_table(ps2), "matrix")
metadata = as(sample_data(ps2), "matrix")
write.table(otu,"Fieldsites_silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"Fieldsites_silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)

#### CLEAR DATA and read back in final cleaned data (no NAs, no chloroplasts, no mitochondria, no eukaryotes)
otu <- read.table("Fieldsites_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Fieldsites_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Fieldsites_metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
#11332 taxa and 65 samples
get_taxa_unique(ps, "Family") #593
get_taxa_unique(ps, "Order") #336
get_taxa_unique(ps, "Kingdom") #2

# remove control samples for plotting, remaining samples = 62
ps = subset_samples(ps, Coral != "control")
ps

# plot number of observed ASVs in coral samples
plot_richness(ps,x="Condition",color="Coral",measures=c("Observed"))

# look at data and chose filtering method for very low abundance ASVs
ntaxa(ps) #11332
ps5<-filter_taxa(ps, function(x) mean(x) >5, TRUE)
ntaxa(ps5) #693
get_taxa_unique(ps, "Genus") # 1254
get_taxa_unique(ps5, "Genus") #279


# filtered ASVs with very low abundance with phyloseq, now export otu and taxa tables from phyloseq for codaseq
otu = as(otu_table(ps5), "matrix")
taxon = as(tax_table(ps5), "matrix")
metadata = as(sample_data(ps5), "matrix")
write.table(otu,"Fieldsites_ps5_silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"Fieldsites_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)
write.table(metadata,"Fieldsites_ps5_silva_metadata.txt",sep="\t",col.names=NA)

######### Perform center-log-ratio transformation on ASVs and calculate Aitchison Distance and principal components
# READ IN OTU data that has been filtered for very low abundance sequences; do not clear data here. Keep phyloseq object ps5 for anosim/permanova
otu <- read.table("Fieldsites_ps5_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Fieldsites_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Fieldsites_ps5_silva_metadata.txt",sep="\t",header=T,row.names=1)
genus<-as.character(taxon$Genus)

# First, replace 0 values with an estimate (because normalization is taking log, can't have 0)
# Also transposing here, need samples as rows
d.czm <- cmultRepl(t(otu), method="CZM", label=0)
# Perform the center-log-ratio (CLR) transformation 
d.clr <- codaSeq.clr(d.czm)
# transpose matrix of CLR transformed data for ordination and dendrogram
E.clr <- t(d.clr)
# plot compositional PCA biplot (perform a singular value decomposition)
d.pcx <- prcomp(E.clr)
# calculate percent variance explained for the axis labels
pc1 <- round(d.pcx$sdev[1]^2/sum(d.pcx$sdev^2),2)
pc2 <- round(d.pcx$sdev[2]^2/sum(d.pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")
biplot(d.pcx, cex=c(0.6,0.4), var.axes=F,scale=1, xlab=xlab, ylab=ylab, ylabs=genus)
summary(d.pcx)
str(d.pcx)
screeplot(d.pcx)

######### FIGURE 2
# replot PCA with ggplot2 (showing samples only)
df_out <- as.data.frame(d.pcx$x)
theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
cols<-c("Dichocoenia stokesii"="#009E73","Diploria labyrinthiformis"="#e79f00","Montastraea cavernosa"="#56B4E9","Orbicella faveolata"="#0072B2")
pdf("Figure2_Fieldsites_PCA_condition_coral.pdf",width=8.5)
p<-ggplot(df_out,aes(x=PC1,y=PC2,fill=samples$Coral,shape=samples$Condition))
p<-p+geom_point(size=3)+theme(axis.title = element_text(size=14))+theme(axis.text=element_text(size=12))+
  theme(legend.title = element_text(size=14))+theme(legend.text = element_text(size=12))+
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,22,24))+
  guides(fill = guide_legend(override.aes=list(shape=21)))
p + labs(x=xlab, y=ylab, fill="Coral", shape="Condition") + coord_fixed()
dev.off()

####### Use phyloseq/vegan to perform ANOSIM/PERMANOVA
# set metadata as factors for anosim
conds<-as.character(samples$Condition)
site<-as.character(samples$Site)
coral<-as.character(samples$Coral)

# anosim between groups using Aitchison distance
dist.clr <- dist(E.clr)
ano.conds <- anosim(dist.clr, conds, permutations=999)
pdf("Fieldsites_ANOSIM_Conds.pdf",width=8.5)
plot(ano.conds)
dev.off()

ano.site <- anosim(dist.clr, site, permutations=999)
pdf("Fieldsites_ANOSIM_Site.pdf",width=8.5)
plot(ano.site)
dev.off()

ano.coral <- anosim(dist.clr, coral, permutations=999)
pdf("Fieldsites_ANOSIM_Coral.pdf",width=8.5)
plot(ano.coral)
dev.off()

# permanova between groups using Aitchison distance
perm<-adonis(dist.clr~conds*site,as(sample_data(ps5),"data.frame"))
print(perm)

perm<-adonis(dist.clr~conds*coral,as(sample_data(ps5),"data.frame"))
print(perm)

########## Figure 3 - bar charts using ps5 (low abundance filtered) with 62 samples (693 taxa)
ps5
ps_ra<-transform_sample_counts(ps5, function(OTU) OTU/sum(OTU))
ps_ra_mcav = subset_samples(ps_ra, Coral == "Montastraea cavernosa")
ps_ra_ofav = subset_samples(ps_ra, Coral == "Orbicella faveolata")
ps_ra_dlab = subset_samples(ps_ra, Coral == "Diploria labyrinthiformis")
ps_ra_dsto = subset_samples(ps_ra, Coral == "Dichocoenia stokesii")
#figure out how many colors you need
get_taxa_unique(ps_ra, "Order") #99
get_taxa_unique(ps_ra, "Class") #33
#you can make n any number of colors you want; with as much difference between the colors as possible (distinct colors)
n <- 99
palette <- distinctColorPalette(n)
#you can rerun the previous line to get a new selection of colors

sample_data(ps_ra_mcav)$Near_Far<-factor(sample_data(ps_ra_mcav)$Near_Far,levels=c("undiseased","far","near","lesion"))
p1=plot_bar(ps_ra_mcav, fill="Class")+
  geom_bar(aes(fill=Class), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  facet_grid(.~Near_Far,scales="free",space="free")+
  ggtitle("Montastraea cavernosa")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank())

sample_data(ps_ra_ofav)$Near_Far<-factor(sample_data(ps_ra_ofav)$Near_Far,levels=c("undiseased","far","near","lesion"))
p2=plot_bar(ps_ra_ofav, fill="Class")+
  geom_bar(aes(fill=Class), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  facet_grid(.~Near_Far,scales="free",space="free")+
  ggtitle("Orbicella faveolata")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position="right")+
  theme(axis.title.x = element_blank())

sample_data(ps_ra_dlab)$Near_Far<-factor(sample_data(ps_ra_dlab)$Near_Far,levels=c("undis.","far","near","lesion"))
p3=plot_bar(ps_ra_dlab, fill="Class")+
  geom_bar(aes(fill=Class), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  facet_grid(.~Near_Far,scales="free",space="free")+
  ggtitle("Diploria labyrinthiformis")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())

sample_data(ps_ra_dsto)$Near_Far<-factor(sample_data(ps_ra_dsto)$Near_Far,levels=c("undiseased","far","near","lesion"))
p4=plot_bar(ps_ra_dsto, fill="Class")+
  geom_bar(aes(fill=Class), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  facet_grid(.~Near_Far,scales="free",space="free")+
  ggtitle("Dichocoenia stokesii")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank())

# adjust width and height until it looks right for double columns
pdf("Figure3_Fieldsites_4species_BarCharts_Class.pdf",width=16, height=10)
plot_grid(p1,p2,p3,p4,labels=c("A","B","C","D"), ncol=2, nrow=2)
dev.off()

########## Figure 4 - Beta Diversity Dispersions
#calculate multivariate dispersions based on condition
mod <-betadisper(dist.clr, conds)
#one way anova
anova(mod)
#boxplots
plot(mod)
boxplot(mod)

## Compute mean distance to centroid per group
#this just prints values on the console
tapply(mod$distances, conds, mean)
## Same, but variance instead
tapply(mod$distances, conds, var)

#Get the distances to centroid from the model
mod$distances
dis <- mod$distances
#melt
dis.melt <- melt(dis)
#move rownames to columns so we can merge the dispersion values and metadata
dis.melt$Sample <- rownames(dis.melt)
samples$Sample <- rownames(samples)
#merge metadata and dispersion 
dis.treat <- merge(samples, dis.melt)
#rename column
colnames(dis.treat)[7] <- "distance"

#run linear model to test significance
distlm <-lm(distance~Coral*Condition, data=dis.treat)
summary(distlm)
anova(distlm)

# plot average dispersion by group, with all points shown
p1<-ggplot(dis.treat, aes(x=Condition,y=distance))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Coral),size=3)+
  scale_color_manual(values=c("#009E73","#e79f00","#56B4E9","#0072B2"))+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(face="italic"))+
  theme(text=element_text(size=16))+
  ylab("Distance to Centroid")
p2<-ggplot(dis.treat, aes(x=Coral,y=distance))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Condition),size=3)+
  scale_color_manual(values=c("#D55E00","#999999","#000000"))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_text(face="italic"))+
  theme(text=element_text(size=14))+
  ylab("Distance to Centroid")
pdf("Figure4_DistanceToCentroid_2panel.pdf",width=11,height=11)
plot_grid(p1,p2,labels=c("A","B"), ncol=1, nrow=2)
dev.off()

############ DESeq2 analysis
# check to see that you are using the full dataset (ps), not the low abundance ASV filtered dataset (ps5)
ps
# you should have 11332 taxa and 62 samples

############ DESeq2 analysis for Mcav, using lesion versus far from lesion samples
ps.mcav = subset_samples(ps, Coral=="Montastraea cavernosa")
ps.mcav
ps.mcav.notnear = subset_samples(ps.mcav, Near_Far != "near" & Near_Far != "undis." & Near_Far != "undiseased")
ps.mcav.notnear
# Define the order of the conditions for testing
# In this order, the positive fold change values are what increased in the disease lesions (DD)
sample_data(ps.mcav.notnear)$Condition<-factor(sample_data(ps.mcav.notnear)$Condition,levels=c("DH","DD"))
head(sample_data(ps.mcav.notnear)$Condition, 10)
# DESEQ2 analysis on MCAV samples
dds = phyloseq_to_deseq2(ps.mcav.notnear, ~ Condition)
dds
#filter rows with very few counts
dds <- dds[ rowSums(counts(dds)) > 5, ]
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds,test="Wald", fitType="parametric")
res = results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
#save table of results
sig = as(sigtab, "matrix")
write.table(sig,"DESeq2_results_Mcav_notnear.txt",sep="\t",col.names=NA)

############ DESeq2 analysis for Dlab, using lesion versus far from lesion samples
ps.dlab = subset_samples(ps, Coral=="Diploria labyrinthiformis")
ps.dlab.notnear = subset_samples(ps.dlab, Near_Far != "near" & Near_Far != "undis." & Near_Far != "undiseased")
ps.dlab.notnear
# Define the order of the conditions for testing
# In this order, the positive fold change values are what increased in the disease lesions (DD)
sample_data(ps.dlab.notnear)$Condition<-factor(sample_data(ps.dlab.notnear)$Condition,levels=c("DH","DD"))
head(sample_data(ps.dlab.notnear)$Condition, 10)
# DESEQ2 analysis on DLAB samples
dds = phyloseq_to_deseq2(ps.dlab.notnear, ~ Condition)
dds
#filter rows with very few counts
dds <- dds[ rowSums(counts(dds)) > 5, ]
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds,test="Wald", fitType="parametric")
res = results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
#save table of results
sig = as(sigtab, "matrix")
write.table(sig,"DESeq2_results_Dlab_notnear.txt",sep="\t",col.names=NA)

############ DESeq2 analysis for Dsto, using lesion versus far from lesion samples
ps.dsto = subset_samples(ps, Coral=="Dichocoenia stokesii")
ps.dsto.notnear = subset_samples(ps.dsto, Near_Far != "near" & Near_Far != "undis." & Near_Far != "undiseased")
ps.dsto.notnear
# Define the order of the conditions for testing
# In this order, the positive fold change values are what increased in the disease lesions (DD)
sample_data(ps.dsto.notnear)$Condition<-factor(sample_data(ps.dsto.notnear)$Condition,levels=c("DH","DD"))
head(sample_data(ps.dsto.notnear)$Condition, 10)
# DESEQ2 analysis on DSTO samples
dds = phyloseq_to_deseq2(ps.dsto.notnear, ~ Condition)
dds
#filter rows with very few counts
dds <- dds[ rowSums(counts(dds)) > 5, ]
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds,test="Wald", fitType="parametric")
res = results(dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)
#save table of results
sig = as(sigtab, "matrix")
write.table(sig,"DESeq2_results_Dsto_notnear.txt",sep="\t",col.names=NA)

######## COMPARE 3 SPECIES (bring them together so that same colors are used for all 3 plots)
dlab<-read.table("DESeq2_results_Dlab_notnear.txt",sep="\t",header=TRUE, row.names=1)
dlab$host<-"dlab"
mcav<-read.table("DESeq2_results_Mcav_notnear.txt",sep="\t",header=TRUE, row.names=1)
mcav$host<-"mcav"
dsto<-read.table("DESeq2_results_Dsto_notnear.txt",sep="\t",header=TRUE, row.names=1)
dsto$host<-"dsto"
dm<-rbind(dlab,mcav)
dmd<-rbind(dm,dsto)
write.table(dmd,"DESeq2_compare3_notnear.txt", sep="\t",col.names=NA)
shared<-read.table("DESeq2_compare3_notnear.txt",sep="\t",header=TRUE)
# pull the species back apart for plotting of shared DE taxa
mcav<-shared[grepl("mcav", shared$host),]
dlab<-shared[grepl("dlab", shared$host),]
dsto<-shared[grepl("dsto", shared$host),]
phyla<-unique(shared$Phylum)
print(phyla)
# total of 24 phyla in shared dataset; print out palette of 24 colors, then assign colors to alphabetically ordered phyla
n <- 24
palette <- distinctColorPalette(n)
print(palette)
cols <- c("Acidobacteria"="#5CE1CF","Actinobacteria"="#DDE594","Bacteria"="#EC8342","Bacteroidetes"="#868296","Chloroflexi"="#E0DACD",
          "Cyanobacteria"="#DE418D","Deferribacteres"="#76DF8E","Epsilonbacteraeota"="#D18278","Euryarchaeota"="#CD67D9",
          "Fibrobacteres"="#813BE7","Firmicutes"="#E7B3C9","Fusobacteria"="#6DC3DB","Halanaerobiaeota"="#6B9FE0",
          "Kiritimatiellaeota"="#E240E9","LCP-89"="#BAECC0","Lentisphaerae"="#DAB279","Marinimicrobia_(SAR406_clade)"="#C496E3",
          "Nanoarchaeaeota"="#758F61","Patescibacteria"="#7567D0","Planctomycetes"="#B2E4E0","Proteobacteria"="#000066",
          "Spirochaetes"="#E0D94E","Thaumarchaeota"="#82E451","Verrucomicrobia"="#C9C5E4")

pdf("FigureS4_DEtaxa_Mcav.pdf",width=8.5,height=11)
p1 <- ggplot(mcav,aes(log2FoldChange,Genus,color=Phylum))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE),height=0.2)+
  scale_color_manual(values=cols)+
  ggtitle("Montastraea cavernosa")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.y=element_blank())+
  theme(axis.title=element_text(size=12))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12))+
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right")
p1
dev.off()

pdf("FigureS5_DEtaxa_Dlab.pdf",width=8.5,height=11)
p2 <- ggplot(dlab,aes(log2FoldChange,Genus,color=Phylum))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE),height=0.2)+
  scale_color_manual(values=cols)+
  ggtitle("Diploria labyrinthiformis")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.y=element_blank())+
  theme(axis.title=element_text(size=12))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12))+
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right")
p2
dev.off()

pdf("FigureS6_DEtaxa_Dsto.pdf",width=8.5,height=15)
p3 <- ggplot(dsto,aes(log2FoldChange,Genus,color=Phylum))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE),height=0.2)+
  scale_color_manual(values=cols)+
  ggtitle("Dichocoenia stokesii")+
  theme(plot.title = element_text(face="italic"))+
  theme(axis.title.y=element_blank())+
  theme(axis.title=element_text(size=12))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=12))+
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right")
p3
dev.off()

############## Figure 5 RELATIVE ABUNDANCE OF DIFFERENTIALLY ABUNDANT ASVs 
otu <- read.table("Fieldsites_ps5_silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("Fieldsites_ps5_silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("Fieldsites_ps5_silva_metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps5 <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps5
ps_ra<-transform_sample_counts(ps5, function(OTU) OTU/sum(OTU))
ps_ra
ps5.subset <- subset_taxa(ps_ra, rownames(tax_table(ps_ra)) %in% c("TACGGAGGGTCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTCCGTAGGCGGGGTTTTAAGTCAGTGGTGAAAGCCTACAGCTCAACTGTAGAACTGCCATTGAAACTGGAACTCTTGAATGTGATTGAGGTAGGCGGAATATGTCATGTAGCGGTGAAATGCTTAGATATGACATAGAACACCGATAGCGAAGGCAGCTTACCAAGTCATTATTGACGCTGATGGACGAAAGCGTGGGGAGCGAACAGG", "TACGTAGGGGGCAAGCGTTATCCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGTTTTTCAAGTCAGAAGTGAAAGGCTATGGCTCAACCATAGTAAGCTTTTGAAACTGTTAAACTTGAGTGCAGGAGAGGAAAGTGGAATTCCTAGTGTAGAGGTGAAATTCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACTGTAACTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGG","TACGGAGGGGGTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATTAGTCAGTCAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTGATACTGCTAGTCTTGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGTGTGGGGAGCAAACAGG","TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTACGCAGGCGGTTAGTTAAGTCAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTTGAAACTGGCTAACTAGAGTGCGACAGAGGGTGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCTGAAGGAATACCGATGGCGAAGGCAGCCACCTGGGTCGACACTGACGCTCATGTACGAAAGCGTGGGTAGCAAACAGG","TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCTGGGCTCAACCCGGGAAGGTCATTTGAAACTGGCAAGCTAGAGTACTGTAGAGGGGGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCTGAAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAGATACTGACACTCAGATGCGAAAGCGTGGGGAGCAAACAGG"))
ps5.subset
otu = as(otu_table(ps5.subset), "matrix")
taxon = as(tax_table(ps5.subset), "matrix")
metadata = as(sample_data(ps5.subset), "matrix")
colnames(otu)[colnames(otu)=="TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCTGGGCTCAACCCGGGAAGGTCATTTGAAACTGGCAAGCTAGAGTACTGTAGAGGGGGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCTGAAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAGATACTGACACTCAGATGCGAAAGCGTGGGGAGCAAACAGG"] <- "Cryomorphaceae"
colnames(otu)[colnames(otu)=="TACGTAGGGGGCAAGCGTTATCCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGTTTTTCAAGTCAGAAGTGAAAGGCTATGGCTCAACCATAGTAAGCTTTTGAAACTGTTAAACTTGAGTGCAGGAGAGGAAAGTGGAATTCCTAGTGTAGAGGTGAAATTCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACTGTAACTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGG"] <- "Fusibacter"
colnames(otu)[colnames(otu)=="TACGGAGGGTCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTCCGTAGGCGGGGTTTTAAGTCAGTGGTGAAAGCCTACAGCTCAACTGTAGAACTGCCATTGAAACTGGAACTCTTGAATGTGATTGAGGTAGGCGGAATATGTCATGTAGCGGTGAAATGCTTAGATATGACATAGAACACCGATAGCGAAGGCAGCTTACCAAGTCATTATTGACGCTGATGGACGAAAGCGTGGGGAGCGAACAGG"] <- "Planktotalea"
colnames(otu)[colnames(otu)=="TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTACGCAGGCGGTTAGTTAAGTCAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTTGAAACTGGCTAACTAGAGTGCGACAGAGGGTGGTAGAATTTCAGGTGTAGCGGTGAAATGCGTAGAGATCTGAAGGAATACCGATGGCGAAGGCAGCCACCTGGGTCGACACTGACGCTCATGTACGAAAGCGTGGGTAGCAAACAGG"] <- "Algicola"
colnames(otu)[colnames(otu)=="TACGGAGGGGGTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATTAGTCAGTCAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTGATACTGCTAGTCTTGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGTGTGGGGAGCAAACAGG"] <- "Vibrio"
otu<-as.data.frame(otu)
otu<-rownames_to_column(otu,var="Sample")
metadata<-as.data.frame(metadata)
metadata<-rownames_to_column(metadata,var="Sample")
otu.meta<-merge(metadata,otu,"Sample")
otu_long<-melt(otu.meta,id.vars=c("Sample","Coral","Colony","Condition","Near_Far","Source","Site"),variable.name="ASV",value.name="Proportion")
otu_long$Coral<-factor(otu_long$Coral,levels=c("Montastraea cavernosa","Orbicella faveolata","Diploria labyrinthiformis","Dichocoenia stokesii"))
asv1<-otu_long[grepl("Cryomorphaceae", otu_long$ASV),]
asv2<-otu_long[grepl("Fusibacter", otu_long$ASV),]
asv3<-otu_long[grepl("Planktotalea", otu_long$ASV),]
asv4<-otu_long[grepl("Algicola", otu_long$ASV),]
asv5<-otu_long[grepl("Vibrio", otu_long$ASV),]

p1<-ggplot(asv1, aes(x=Source,y=Proportion))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Source),size=3)+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=14))+
  facet_grid(.~Coral)+
  theme(strip.text.x=element_text(face="italic",size=10))+
  ylab("Relative Abundance")+
  ggtitle("Cryomorphaceae")

p2<-ggplot(asv2, aes(x=Source,y=Proportion))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Source),size=3)+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=14))+
  facet_grid(.~Coral)+
  theme(strip.text.x=element_text(face="italic",size=10))+
  ylab("Relative Abundance")+
  ggtitle("Fusibacter")+
  theme(plot.title = element_text(face="italic"))

p3<-ggplot(asv3, aes(x=Source,y=Proportion))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Source),size=3)+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=14))+
  facet_grid(.~Coral)+
  theme(strip.text.x=element_text(face="italic",size=10))+
  ylab("Relative Abundance")+
  ggtitle("Planktotalea")+
  theme(plot.title = element_text(face="italic"))

p4<-ggplot(asv4, aes(x=Source,y=Proportion))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Source),size=3)+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=14))+
  facet_grid(.~Coral)+
  theme(strip.text.x=element_text(face="italic",size=10))+
  ylab("Relative Abundance")+
  ggtitle("Algicola")+
  theme(plot.title = element_text(face="italic"))
  
p5<-ggplot(asv5, aes(x=Source,y=Proportion))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0),aes(color=Source),size=3)+
  theme(axis.title.x=element_blank())+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=14))+
  facet_grid(.~Coral)+
  theme(strip.text.x=element_text(face="italic",size=10))+
  ylab("Relative Abundance")+
  ggtitle("Vibrio")+
  theme(plot.title = element_text(face="italic"))

pdf("Figure5_RelAbund_5ASVs.pdf",width=8.5,height=11)
plot_grid(p1,p2,p3,p4,p5, labels=c("A","B","C","D","E"), ncol=1, nrow=5)
dev.off()

                                    
###### Figures S1, S2, S3 - Bar charts with one coral species at a time, finer resolution than Class
get_taxa_unique(ps_ra_mcav, "Genus") #279
# get rid of ASVs with no counts in mcav
ps_ra_mcav <- prune_taxa(taxa_sums(ps_ra_mcav) != 0, ps_ra_mcav)
get_taxa_unique(ps_ra_mcav, "Genus") #253
get_taxa_unique(ps_ra_mcav, "Family") #155

n <- 155
palette <- distinctColorPalette(n)
sample_data(ps_ra_mcav)$Near_Far<-factor(sample_data(ps_ra_mcav)$Near_Far,levels=c("undiseased","far","near","lesion"))
pdf("FigureS1_Mcav_Families.pdf",width=16, height=10)
p1=plot_bar(ps_ra_mcav, fill="Family")+
  geom_bar(aes(fill=Family), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  facet_grid(.~Near_Far,scales="free",space="free")+
  ggtitle("Montastraea cavernosa")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position = "bottom")+
  theme(axis.title.x = element_blank())
p1
dev.off()

get_taxa_unique(ps_ra_dlab, "Genus") #279
# get rid of ASVs with no counts in dlab
ps_ra_dlab <- prune_taxa(taxa_sums(ps_ra_dlab) != 0, ps_ra_dlab)
get_taxa_unique(ps_ra_dlab, "Genus") #216
get_taxa_unique(ps_ra_dlab, "Family") #139
n <- 139
palette <- distinctColorPalette(n)
sample_data(ps_ra_dlab)$Near_Far<-factor(sample_data(ps_ra_dlab)$Near_Far,levels=c("undis.","far","near","lesion"))
pdf("FigureS2_Dlab_Families.pdf",width=16, height=10)
p3=plot_bar(ps_ra_dlab, fill="Family")+
  geom_bar(aes(fill=Family), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  facet_grid(.~Near_Far,scales="free",space="free")+
  ggtitle("Diploria labyrinthiformis")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())
p3
dev.off()

get_taxa_unique(ps_ra_dsto, "Genus") #279
# get rid of ASVs with no counts in dsto
ps_ra_dsto <- prune_taxa(taxa_sums(ps_ra_dsto) != 0, ps_ra_dsto)
get_taxa_unique(ps_ra_dsto, "Genus") #224
get_taxa_unique(ps_ra_dsto, "Family") #139
n <- 139
palette <- distinctColorPalette(n)
sample_data(ps_ra_dsto)$Near_Far<-factor(sample_data(ps_ra_dsto)$Near_Far,levels=c("undis.","far","near","lesion"))
pdf("FigureS3_Dsto_Families.pdf",width=16, height=10)
p4=plot_bar(ps_ra_dsto, fill="Family")+
  geom_bar(aes(fill=Family), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+
  theme(axis.text.x=element_text(angle = 90))+scale_fill_manual(values=palette)+
  facet_grid(.~Near_Far,scales="free",space="free")+
  ggtitle("Dichocoenia stokesii")+
  theme(plot.title = element_text(face="italic"))+
  theme(legend.position="bottom")+
  theme(axis.title.x = element_blank())
p4
dev.off()

