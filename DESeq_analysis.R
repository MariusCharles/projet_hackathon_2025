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

