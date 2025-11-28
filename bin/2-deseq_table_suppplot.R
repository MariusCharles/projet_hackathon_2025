#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)

# Récupérer les arguments
args <- commandArgs(trailingOnly = TRUE)
counts_file  <- args[1]
coldata_file <- args[2]

# Lire counts matrix 
cts <- read.table(
  counts_file,
  header       = TRUE,
  comment.char = "#",
  row.names    = 1,
  check.names  = FALSE
)

# === Pré-traitement ===
# Retirer les colonnes d'annot
annot_cols <- c("Chr", "Start", "End", "Strand", "Length")
cts <- cts[, !(colnames(cts) %in% annot_cols), drop = FALSE]
# Clean noms de colonne 
old_names <- colnames(cts)
new_names <- gsub("_trimmed.sorted.bam$", "", old_names)
colnames(cts) <- new_names
# =====================

# Lire le fichier config
coldata <- read.csv(
  coldata_file,
  row.names   = 1,
  check.names = FALSE
)

# On garde que la colonne 'condition'
coldata <- coldata["condition"]
coldata$condition <- factor(coldata$condition)

# === Quality Check ===     # Aucun output = OK. 
common_samples <- intersect(colnames(cts), rownames(coldata))
if (length(common_samples) != ncol(cts)) {
  stop(
    "Not all samples match.\n",
    "Samples in  counts: ", paste(colnames(cts), collapse = ", "), "\n",
    "Samples in coldata: ", paste(rownames(coldata), collapse = ", "), "\n"
  )
}
# =====================

# Même ordre 
cts <- cts[, rownames(coldata), drop = FALSE]


# === DESEQ2 ===
# Build object DESEQ2
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)
# =============


# Save file 
write.csv(as.data.frame(res), file = "deseq2_results.csv", row.names = TRUE)

# Plot supplementary figure
pdf("MA_plot_allgenes.pdf")
plotMA(res, ylim = c(-4,4), alpha = 0.05)
dev.off()


# === Génération du plot PCA ===
transformed_counts = vst(dds, blind = FALSE)

pdf("PCA_repro.pdf")
plotPCA(transformed_counts, intgroup="condition") +
  labs(color = "Sample",
       title = "PCA - Repro") +
  geom_label(
    aes(label = name),
    size = 3,
    label.padding = unit(0.15, "lines"),
    alpha = 0.7,
    nudge_y = 1.1
  ) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(clip = "off")             # labels en dehors visibles
dev.off()
