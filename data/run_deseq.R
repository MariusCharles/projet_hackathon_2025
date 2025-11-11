#!/usr/bin/env Rscript
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
counts_file  <- args[1]   # counts_matrix.txt (de featurecounts)
coldata_file <- args[2]   # config.csv

## lire mat 
cts <- read.table(
  counts_file,
  header       = TRUE,
  comment.char = "#",
  row.names    = 1,
  check.names  = FALSE
)

## On enlève les colonnes de coordonnées (Chr, Start, End, Strand, Length)
cts <- cts[, !(colnames(cts) %in% c("Chr", "Start", "End", "Strand", "Length"))]

## Clean pour enlever suffixe .sorted.bam 
colnames(cts) <- sub("\\.sorted\\.bam$", "", colnames(cts))
colnames(cts) <- sub("\\.bam$", "", colnames(cts))

## Lire le config file
coldata <- read.csv(
  coldata_file,
  row.names   = 1,
  check.names = FALSE
)

# Garder que la colonne sample et la col condition 
coldata <- coldata[, "condition", drop = FALSE]

# S'assurer que les colonnes de counts correspondent aux samples de coldata
common_samples <- intersect(colnames(cts), rownames(coldata))

# même ordre 
cts     <- cts[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

## construire datasetdeseq
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = coldata,
                              design    = ~ condition)

dds <- DESeq(dds)
res <- results(dds) # run deseq

## Save 
write.csv(as.data.frame(res), file = "deseq2_results.csv", row.names = TRUE)

# Premier MAplot test 
pdf("MA_plot.pdf")
plotMA(res, ylim = c(-4,4), alpha = 0.05)
dev.off()



