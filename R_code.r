####################################################################################################################################################
# This R-code related to the analyzes performed in paper entitled "Unraveling the Molecular Heterogeneity in type 2 Diabetes: A Potential            #
#Subtype Discovery Study Followed by Metabolic Modeling                                                                                            #
# Maryam Khoshnejat1, Kaveh Kavousi1, Ali Mohammad Banaei- Moghaddam, Ali Akbar Moosavi-Movahedi                                                   #
# Laboratory of Complex Biological Systems and Bioinformatics (CBB), Department of Bioinformatics, Institute of Biochemistry and Biophysics (IBB), #
# University of Tehran, Tehran, Iran                                                                                                               #
# Corresponding author: Kaveh Kavousi (kkavousi@ut.ac.ir)                                                                                          #
####################################################################################################################################################



################################################################################## Differential gene expression analysis using Deseq2 package in bioconductor R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2") 
library(DESeq2)


#read data file
cts <- as.matrix(read.csv("t2d_ngt_count_data.txt", row.names = 1, header= TRUE, sep="\t"))
coldata <- read.csv("coldata.txt", row.names = 1, sep="\t")
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)


#filter low count read  
keep <- rowSums(counts(dds)>=5) >= 38
dds <- dds[keep,]


dds$condition <- factor(dds$condition, levels = c("NGT","T2D"))

##### 
dds <- DESeq(dds)
res <- results(dds)


# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

#save result
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), file="results.csv")




########################################################################################### Normalizing gene expression
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
n<-counts(dds, normalized=TRUE)
 
 
 
%%%%%% deseq2 normalization with gene length adjustment

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("AnnotationDbi")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicFeatures")

library(GenomicFeatures)
txdb <- makeTxDbFromGFF("gencode.v19.annotation.gtf_withproteinids", format="gtf", circ_seqs=character())
ebg <- exonsBy(txdb, by="gene") 
rowRanges(dds) = GRangesList(ebg)
FPKM = fpkm(dds)
 
###########################################################################################get gene name and GO term from ensmbl id
mygenes <- read.csv("deseq_genes_id.txt", header= TRUE, sep="\t")
genes <- mygenes$id
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","hgnc_id", "go_id", "namespace_1003","name_1006","description"),values=genes,mart= ensembl)

totaltable<-merge(mygenes,G_list,by.x="id",by.y="ensembl_gene_id")
write.csv(as.data.frame(totaltable), file="Go_terms.csv")

############################################################################################### Hierarchical clustering
cts <- t(as.matrix(read.csv("T2d_exp.txt", row.names = 1, header= TRUE, sep="\t")))
install.packages ("hyperSpec")
library ("hyperSpec")
d <- dist(cts, method = 'euclidean')
hc1 <- hclust(d, method = "complete" )

install.packages("dendextend")
install.packages("circlize")
library(dendextend)
library(circlize)

# create a dendrogram
dend <- as.dendrogram(hc1)
# change the dendrogram to have colors in the branches 
dend <- dend %>% 
  color_branches(k=3) %>% 
  color_labels

# the radial plot
par(mar = rep(0,4))
# circlize_dendrogram(dend, dend_track_height = 0.8) 
circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .3) 


################################################################################################## Heatmap
#pheatmap
library(pheatmap)
cts <-t(as.matrix(read.csv("DEG.txt", row.names = 1, header= TRUE, sep="\t")))
pheatmap((cts), cluster_rows=FALSE, cluster_cols=TRUE, scale ='column' , border_color="NA", cellwidth = 10, cellheight =30, col=colors)

