#Library load
install.packages("DESeq2")
library(DESeq2)
library(tximport)
library(readr)
library(rhdf5)
library(rlang)
library(dplyr)
install.packages("tidyverse")
library(tidyverse)
install.packages("biomaRt")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt", force = TRUE)
library(biomaRt)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("psych")
library(psych)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(RColorBrewer)
setwd("~/Desktop/dep_lncrna/")
samples <- read.csv("~/Desktop/dep_lncrna/meta_data.csv")
files <- list.files(path = "~/Desktop/dep_lncrna/New_Kallisto_out/", pattern = "*.tsv$", recursive=TRUE, full.names = TRUE)
names(files) <- samples$id
txi <- tximport(files = files, type = "kallisto", txOut = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ treatment)
ds <- DESeq(ddsTxi, test="Wald")
df <- results(ds, contrast = c("treatment", "COPD_control", "COPD"))
healthy_df <- results(ds, contrast = c("treatment", "Healthy_control", "Healthy"))
#df1 <- data.frame(df) %>% filter(padj < 0.05)
df1 <- add_rownames(data.frame(df), var = "Gene_id")
head(df)
head(df1)
healthy_df1 <- add_rownames(data.frame(healthy_df), var = "Gene_id")

anno_transcript <- read.csv("~/Desktop/dep_lncrna/annotation_file.tsv", sep = "\t", header = F)
#anno_transcript <- read.csv("~/Desktop/dep_lncrna/transcript_annotation.tsv", sep = "\t", header = F)
#tail(anno_transcript)
#colnames(anno_transcript)
head(anno_transcript, 10)
colnames(anno_transcript) <- c("Transcript_id", "Gene_id", "Transcript_class", "Type")
#colnames(anno_transcript) <- c("Transcript_id", "Gene_id", "Transcript_class")
df2 <- merge(df1, anno_transcript, by="Gene_id", all=F)
head(df2)
tail(df2)
mRNA <- filter(df2, Transcript_class=="protein_coding")
lncRNA <- filter(df2, Transcript_class=="lncRNA")
#for healthy
healthy_df2 <- merge(healthy_df1, anno_transcript, by="Gene_id", all=F)
healthy_mRNA <- filter(healthy_df2, Transcript_class=="protein_coding")
healthy_lncRNA <- filter(healthy_df2, Transcript_class=="lncRNA")
#volcano plot for COPD
with(lncRNA, plot(log2FoldChange, -log10(pvalue), pch=20, main="C vs COPD (lncRNA)", xlim=c(-25,24), cex = 1.5, cex.axis=1.5, cex.lab=1.5))
with(subset(lncRNA, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex = 1.5))
with(subset(lncRNA, padj<.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex = 1.5))

with(mRNA, plot(log2FoldChange, -log10(pvalue), pch=20, main="C vs COPD (mRNA)", xlim=c(-25,24), cex = 1.5, cex.axis=1.5, cex.lab=1.5))
with(subset(mRNA, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex = 1.5))
with(subset(mRNA, padj<.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex = 1.5))

#volcano plot for Healthy
with(healthy_lncRNA, plot(log2FoldChange, -log10(pvalue), pch=20, main="C vs Healthy (lncRNA)", xlim=c(-25,24), cex = 1.5, cex.axis=1.5, cex.lab=1.5))
with(subset(healthy_lncRNA, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex = 1.5))
with(subset(healthy_lncRNA, padj<.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex = 1.5))

with(healthy_mRNA, plot(log2FoldChange, -log10(pvalue), pch=20, main="C vs Healthy (mRNA)", xlim=c(-25,24), cex = 1.5, cex.axis=1.5, cex.lab=1.5))
with(subset(healthy_mRNA, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex = 1.5))
with(subset(healthy_mRNA, padj<.05 & log2FoldChange < -1), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex = 1.5))

#Saving log2FC file as csv file for COPD
write.csv(mRNA, file = "~/Desktop/dep_lncrna/DESeq2 output_20211126/CvsCOPD_mRNA.csv")
write.csv(lncRNA, file = "~/Desktop/dep_lncrna/DESeq2 output_20211126/CvsCOPD_lncRNA.csv")

#Saving log2FC file as csv file for Healthy
healthy_lncRNA <- subset(healthy_lncRNA, padj<.05 & abs(log2FoldChange) > 1)
healthy_mRNA <- subset(healthy_mRNA, padj<.05 & abs(log2FoldChange) > 1)

write.csv(healthy_mRNA, file = "C~/Desktop/dep_lncrna/new_CvsHealthy_mRNA.csv")
write.csv(healthy_lncRNA, file = "~/Desktop/dep_lncrna//new_CvsHealthy_lncRNA.csv")



#GSEA analysis
library(fgsea)
#Gene set enrichment analysis of cancer related genes
negative_lncRNA <- read.table("C:/Users/Lin Lab/Desktop/Sabbir-3/Computational_Biology/DEP_Project/negative.lncRNA.glist.xls.txt")
positiveative_lncRNA <- read.table("C:/Users/Lin Lab/Desktop/Sabbir-3/Computational_Biology/DEP_Project/positive.lncRNA.glist.xls.txt")
all_lncRNA <- rbind(negative_lncRNA, positiveative_lncRNA)
colnames(all_lncRNA) <- c("Gene_id")
genesets <- list(Cancer_Related=all_lncRNA$Gene_id)
head(positiveative_lncRNA)
#gsea_analysis <- function()
gseaDat <- lncRNA %>% filter(!is.na(Gene_id)) %>% filter(padj < 0.05)
ranks <- gseaDat$padj
names(ranks) <- gseaDat$Gene_id
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)

gseaDat <- lncRNA %>% filter(!is.na(Gene_id)) %>% filter(padj < 0.05)
ranks <- gseaDat$padj
names(ranks) <- gseaDat$Geneid1
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=4, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)

gseaDat <- CvsOH %>% filter(!is.na(Geneid1)) %>% filter(log2FoldChange>0)
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$Geneid1
ranks <- sort(ranks, decreasing = T)
barplot(sort(ranks, decreasing = T))
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)

gseaDat <- CvsOH %>% filter(!is.na(Geneid1)) %>% filter(log2FoldChange<0)
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$Geneid1
ranks <- sort(ranks, decreasing = T)
barplot(sort(ranks, decreasing = T))
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)

--------------------------------------------------------------------------------------------

  #Analysis of gene regulated by lncRNA in trans
CvsO3_lncRNA <- CvsO3 %>% filter(biotype=="lncRNA") %>% arrange(padj)
CvsO3_lncRNA_20 <- CvsO3_lncRNA[1:25,]
write.csv(CvsO3_lncRNA_20, "./CvsO3_lncRNA_25.csv", quote = F)
CvsOH_lncRNA <- CvsOH %>% filter(biotype=="lncRNA") %>% arrange(padj)
CvsOH_lncRNA_20 <- CvsOH_lncRNA[1:25,]
write.csv(CvsOH_lncRNA_20, "./CvsOH_lncRNA_25.csv", quote = F)

lncRNA_mRNA <- read.csv("./lncRNA_RNA interaction_final.csv", header = TRUE)
lncRNA_mRNA <- lncRNA_mRNA[complete.cases(lncRNA_mRNA), ] %>% rename(GeneSymble=Name)
CvsO3_interaction <- merge(CvsO3, lncRNA_mRNA, by="GeneSymble", all=F) %>% rename(mRNA=Geneid)
CvsO3_interaction <- CvsO3_interaction %>% rename(Geneid=Name_lncRNA)
CvsO3_interaction <- merge(CvsO3_interaction, id_symbol, by="Geneid", all=F)

CvsO3_interaction1 <- CvsO3_interaction %>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
write.csv(CvsO3_interaction1, "./CvsO3_network.csv", quote = F)
CvsO3_interaction1 <- CvsO3_interaction1[,c(1,7)]

CvsOH_interaction <- merge(CvsOH, lncRNA_mRNA, by="GeneSymble", all=F) %>% rename(mRNA=Geneid)
CvsOH_interaction <- CvsOH_interaction %>% rename(Geneid=Name_lncRNA)
CvsOH_interaction <- merge(CvsOH_interaction, id_symbol, by="Geneid", all=F)
CvsOH_interaction1 <- CvsOH_interaction[,c(2,5,9,10,14,18,20,21)]
CvsOH_interaction1 <- CvsOH_interaction1 %>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
write.csv(CvsOH_interaction1, "./CvsOH_network.csv", quote = F)
CvsOH_interaction1 <- CvsOH_interaction1[,c(1,7)]

library(GGally)
library(ggnet)
library(ggnetwork)
library(network)
#library(statnet.common)
library(sna)
library(ggplot2)
setwd("C:/Users/Lin Lab/Desktop/Sabbir-3/Computational_Biology/DEP_Project/")
lncRNA_mRNA <- read.csv("C:/Users/Lin Lab/Desktop/Sabbir-3/Computational_Biology/DEP_Project/100_lncRNA-RNA_interaction.csv", header = TRUE)
lncRNA_mRNA <- lncRNA_mRNA[complete.cases(lncRNA_mRNA), ] # %>% rename(GeneSymble=Name)
DEG <- read.csv("C:/Users/Lin Lab/Desktop/Sabbir-3/Computational_Biology/DEP_Project/DESeq2 output_051421/new_CvsHealthy_mRNA.csv", header = TRUE)
head(lncRNA_mRNA)
head (distinct(DEG))



levs <- unique(unlist(lncRNA_mRNA, use.names = FALSE))
net <- table(lapply(lncRNA_mRNA, factor, levs))
label <- colnames(net)
net = network(net, directed = FALSE)
attribute <- as.data.frame(label)
colnames(attribute) <- c("GeneSymble")
OH <- inner_join(attribute, CvsOH)
OH$status <- factor(ifelse(OH$log2FoldChange >0, "Up", "Down")) 
network.vertex.names(net) = label
net %v% "Status" <- as.character(OH$status)
color = c(rep(c("mRNA"), 203),rep(c("lncRNA"), 20))
net %v% "Type" <- color
ggnet2(net, node.size = 4, label = TRUE, label.size=2.5, 
       color = "Type", shape = "Status", palette = "Set1")


levs <- unique(unlist(CvsOH_interaction1, use.names = FALSE))
net <- table(lapply(CvsOH_interaction1, factor, levs))
label <- colnames(net)
net = network(net, directed = FALSE)
attribute <- as.data.frame(label)
colnames(attribute) <- c("GeneSymble")
OH <- inner_join(attribute, CvsOH)
OH$status <- factor(ifelse(OH$log2FoldChange >0, "Up", "Down")) 
network.vertex.names(net) = label
net %v% "Status" <- as.character(OH$status)
color = c(rep(c("mRNA"), 201),rep(c("lncRNA"), 20))
net %v% "Type" <- color
ggnet2(net, node.size = 4, label = TRUE, label.size=2.5, 
       color = "Type", shape = "Status", palette = "Set1")

CvsO3_lncRNA_20 <- CvsO3_lncRNA_20 %>% rename(GeneSymble.y=GeneSymble)
CvsO3_interaction <-  CvsO3_interaction%>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
SI_Table4 <- merge(CvsO3_interaction, CvsO3_lncRNA_20, by="GeneSymble.y", all=F)
write.csv(SI_Table4, "SI_Table4.csv")
CvsOH_lncRNA_20 <- CvsOH_lncRNA_20 %>% rename(GeneSymble.y=GeneSymble)
CvsOH_interaction <-  CvsOH_interaction%>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
SI_Table4 <- merge(CvsO3_interaction, CvsO3_lncRNA_20, by="GeneSymble.y", all=T)
write.csv(SI_Table4, "SI_Table5.csv")


