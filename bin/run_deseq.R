#!/usr/bin/env Rscript

library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
counts_file  <- args[1]   # counts_matrix.txt
coldata_file <- args[2]   # config.csv

# Lire counts matrix 
cts <- read.table(
  counts_file,
  header       = TRUE,
  comment.char = "#",
  row.names    = 1,
  check.names  = FALSE
)

# Retirer les colonnes d'annot
annot_cols <- c("Chr", "Start", "End", "Strand", "Length")
cts <- cts[, !(colnames(cts) %in% annot_cols), drop = FALSE]

# Clean noms de colonne
colnames(cts) <- basename(colnames(cts))

colnames(cts) <- sub("\\.sorted\\.bam$", "", colnames(cts))
colnames(cts) <- sub("\\.bam$", "", colnames(cts))

colnames(cts) <- sub("_trimmed$", "", colnames(cts))

# Lire le fichier config
coldata <- read.csv(
  coldata_file,
  check.names = FALSE
)

# sample en rownames
rownames(coldata) <- coldata$sample

# on garde que la colonne condition
coldata <- coldata["condition"]


# petit checkup (debuggage)
common_samples <- intersect(colnames(cts), rownames(coldata))

if (length(common_samples) == 0) {
  stop(
    "Aucun échantillon en commun entre counts_matrix et config.csv.\n",
    "Noms dans counts: ", paste(colnames(cts), collapse = ", "), "\n",
    "Noms dans coldata: ", paste(rownames(coldata), collapse = ", "), "\n"
  )
}

cts     <- cts[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

# build objet deseq2
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

# Save file + plot figure supp
write.csv(as.data.frame(res), file = "deseq2_results.csv", row.names = TRUE)

pdf("MA_plot.pdf")
plotMA(res, ylim = c(-4,4), alpha = 0.05)
dev.off()
