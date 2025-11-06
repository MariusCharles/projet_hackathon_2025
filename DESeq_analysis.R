library("DESeq")

cts <- as.matrix(read.csv("counts_matrix.csv", row.names = 1))
coldata <- read.csv("coldata.csv", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = cts,  # matrice cts directement issue de FeatureCounts
                              colData = coldata,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds)


dds <- DESeq(dds)

# résultats extraits
res <- results(dds)

# save table de résultats
write.csv(as.data.frame(res), file = "deseq2_results.csv")

# open table de résultat
#res=read.csv("deseq2_results.csv")

# Check
print(head(res))

# Keep only <0.05 pval genes
de_genes <- rownames(res)[which(res$padj < 0.05)]

pdf("MA_plot.pdf")
plotMA(de_genes, ylim=c(-4,4))
dev.off()
