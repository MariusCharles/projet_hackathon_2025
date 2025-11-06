library("DESeq2")

# Restructurer matrice de comptes
# Créer Coldata persister/control 

cts <- as.matrix(read.csv("counts_matrix.csv", row.names = 1))
coldata <- read.csv("coldata.csv", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = cts,  # matrice cts directement issue de FeatureCounts
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

# save table de résultats
write.csv(as.data.frame(res), file = "deseq2_results.csv")

# open table de résultat
#res=read.csv("deseq2_results.csv")

# Check
print(head(res))

# Plot SuppFig3
pdf("MA_plot.pdf")
plotMA(res, ylim=c(-4,4),alpha=0.05)
dev.off()

# Plot Fig3
