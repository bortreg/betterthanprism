## Import Normalized mRNA counts from a Nanostring experiment for analysis in Bioconductor DESeq2 package
 # This analysis is desingned to compare Fold Changes across 2 different groups 
 # (eg. Group A and B both receive treatment X, compare responses to X in A vs B)

library(data.table)


#Import and transpose data
dfA <- read.csv("~/Documents/R/Nanostring_DataFrames/HInormDESeq.csv")
dfB <- read.csv("~/Documents/R/Nanostring_DataFrames/LOnormDESeq.csv")
row.names(dfA) <- dfA$sample.ID
row.names(dfB) <- dfB$sample.ID
dfA$sample.ID <- NULL
dfB$sample.ID <- NULL
dfA_trans <- transpose(dfA)
dfB_trans <- transpose(dfB)
rownames(dfA_trans) <- colnames(dfA)
colnames(dfA_trans) <- rownames(dfA)
rownames(dfB_trans) <- colnames(dfB)
colnames(dfB_trans) <- rownames(dfB)
dfA <- dfA_trans
dfB <- dfB_trans

#Run example to test DESeq using matirx of data with existing counts
library(pasilla)
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]