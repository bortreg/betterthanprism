## Import Normalized mRNA counts from a Nanostring experiment for analysis in Bioconductor DESeq2 package
 # This analysis is designed to compare Fold Changes across 2 different groups 
 # (eg. Group A and B both receive treatment X, compare responses to X in A vs B)
 # Followed Description at https://support.bioconductor.org/p/69705/
 # Tutorial for DESeq2 can be found at http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library(data.table)
library(ggplot2)

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
countData <- cbind(dfA, dfB)
countData <- round(countData)

#Create coldata object for DESeq2
colData <- data.frame(
  row.names = colnames(countData),
  "position" = rep(c("A","B"), each=(ncol(countData)/2)),
  "condition" = rep(c("untreated", "treated"), times=(ncol(countData)/2))
)

#check if colData matches countData
all(rownames(colData) %in% colnames(countData))
all(rownames(colData) == colnames(countData))
ncol(countData) == nrow(colData)

#Initial DE setup in DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(countData),
                              colData = colData,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10 #Prefilter reads
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("untreated","treated")) #Set factor levels
dds$position <- factor(dds$position, levels = c("A","B")) 

#Run Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds, name = "condition_treated_vs_untreated")
res <- results(dds, contrast = c("condition","treated","untreated"))
res

write.csv(res, "Desktop/res.csv")

#change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")

#reorder results by p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),] #reorder by p value
sum(res$padj <0.1, na.rm = TRUE) #

#Plot Data
vplot <- res[order(res$padj, decreasing = FALSE),]
pdf("~/Desktop/test1.pdf", width=5, height=5)
ggplot(res, aes(x=log2FoldChange, y=-log2(pvalue), color = FDR)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("#57C84D","#83D475", "#ABE098", "#C5E8B7")) +
  theme_bw() +
  geom_hline(yintercept = -log2(0.0320), linetype = "solid", size =0.75, color = "darkgray") +
  geom_hline(yintercept = -log2(0.00375), linetype = "dashed", size =0.75, color = "darkgray") +
  geom_hline(yintercept = -log2(0.00127), linetype = "dotdash", size =0.75, color = "darkgray") +
  geom_text_repel(
    data = dfA[1:100,],
    aes(label = Gene),
    color = "black",
    size = 2,
    segment.size = 0.1,
    box.padding = unit(0.2, "lines"),
    point.padding = unit(0.25, "lines")
  )
dev.off()


##Run example to test DESeq using matirx of data with existing counts
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

#Adjust coldata to fit cts column names and order
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
cts <- cts[,rownames(coldata)]
all(rownames(coldata) == colnames(cts))

#Run DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

#Add additional feature data
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

#Set factor levels for comparisons of conditions.
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

#Differential Expresion Analysis!
dds <- DESeq(dds)
res <- results(dds)

#Specify coefficient or contrast
res <- results(dds, name="condition_treated_vs_untreated")
res <- results(dds, contrast=c("condition","treated","untreated"))

#Rank by Log Fold Change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")





