#Differential gene expression for bulkRNAseq data

library(DESeq2)
library(gplots)
library(plyr)
library(tidyverse)

# Read in the metadata
metadata<-read.csv("bulk_metadata.csv", header=TRUE, row.names = 1)

#Read in the genecounts
genecounts <- read.csv("bulk_genecounts.csv", header=TRUE, check.names = FALSE, row.names = 1) 

#Make gene counts integers and reorder columns to match metadata
genecountspc <- apply(genecounts, c(1,2), function(x) { (as.integer(x))}) 
genecountspc <- genecountspc[, match(rownames(metadata), colnames(genecountspc))]


#FILTERING
#Capture genes present in at least 30% of samples.
mabsent <- c()

for (i in 1:ncol(genecountspc)) {
        ## Zero spot =absent
        absent <- as.numeric(genecountspc[,i]) < 1
        ## Add the sample to an "absent" matrix with all other samples
        mabsent <- cbind(mabsent, absent)
}

## Determine which probes are present
mpresent <- mabsent==0

## Find which probes are present in at least 30% of the samples
mfilter<- apply(mpresent,1,sum) > ncol(genecountspc)*.3
mpresent30<-genecountspc[mfilter=="TRUE",]
dim(mpresent30)


###DESeq
dds<-DESeqDataSetFromMatrix(countData= mpresent30, colData= metadata, design=~Groups)

dds<-DESeq(dds)
resultsNames(dds)
res<-results(dds)
xres<-res[order(res$padj),]

#Pick a contrast
res<-results(dds, contrast=c("Groups", "VAP_early", "No_VAP_early"))

#Reorder results by log2foldchange
res <- res[order(res$log2FoldChange ), ]

#Add gene symbols
res <- data.frame("gene_symbol" = genecounts[rownames(res), "gene_symbol"], res)

#pull out just FDR<0.05
resSig.05 <- res[which(res$padj < 0.05), ]


#HEATMAPS

#regularized log transformation
rld <- rlog(dds)
rldassay<-assay(rld)

#Get top 50 genes by padj
sigind <- slice_min(res, n=50, order_by = padj)
sigind<-match(rownames(sigind), rownames(rldassay))
set<-rldassay[sigind,]
set <- data.frame("gene_symbol" = genecounts[rownames(set), "gene_symbol"], set, check.names = FALSE)
rownames(set) <- set[, 1]
set <- set[, -1]

#t(set) transposes to allow for centering and scaling
sett<-t(set)
sets<-scale(sett, center=T, scale=T)                   
setst<-t(sets)
rm(sett)   

dev.off()
library("RColorBrewer")
library("gplots")

#Heatmap colors
brewer=rgb(colorRamp(c ("deepskyblue3", "#F9F4EC", "deeppink4"), space="rgb", interpolate="linear")(0:255/255), maxColorValue=255)

library(pheatmap)
metadata_filtered <- metadata[metadata$Groups == "VAP_early"| metadata$Groups == "No_VAP_early", ]
annotation <- metadata_filtered[, 4, drop = FALSE]
colors <- list(Groups = c(VAP_early="deepskyblue1", No_VAP_early="palegreen3"))
setst <- setst[, colnames(setst) %in% rownames(metadata_filtered)]

heatmap <- pheatmap(setst, color = brewer, show_rownames = T, show_colnames = F,
                    scale = "row", border_color = NA, fontsize_col = 10,
                    fontsize_row = 8, cluster_rows = T, cluster_col = T, cutree_row = 1, cutree_col = 1, treeheight_row = 20, treeheight_col = 20,
                    annotation_legend = TRUE, annotation_col = annotation, annotation_colors = colors[1],
                    clustering_distance_cols = "euclidean", clustering_distance_rows = "correlation",
                    cellwidth = 8, cellheight = 8)
