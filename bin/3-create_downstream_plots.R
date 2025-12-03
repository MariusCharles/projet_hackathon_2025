#!/usr/bin/env Rscript
library(ggplot2)
library(ggrepel)

# ==== Arguments Nextflow ====
args <- commandArgs(trailingOnly = TRUE)
deseq_results_path  <- args[1]
genes_pathways_path <- args[2]
mapping_path        <- args[3]

# Lire résultats deseq
res <- read.csv(deseq_results_path,row.names =1)

# Lire table kegg 
kegg_all <- read.table(
  genes_pathways_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Lire table de mapping
mapping <- read.table(
  mapping_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = F
)

# On ajoute 'pth' à la main car il n'est pas présent dans la table d'origine
mapping[mapping$product == "peptidyl-tRNA hydrolase", ]$symbol="pth"


# === Préparation du dataframe pour le plot : ===
# Compute le log2base mean pour les abscisses (si diff de 0)
res$log2BaseMean <- ifelse(res$baseMean > 0, log2(res$baseMean), NA)
# colonne significatif / non significatif 
res$is_sig <- !is.na(res$padj) & res$padj < 0.05
res$signif <- ifelse(res$is_sig, "Significant", "Non-Significant")
# ===============================================


# == On merge les dataframes : ==
# On concatène les colonnes du mapping sur les résultats DESeq
res$locus_tag <- rownames(res)
res_annot <- merge(res, mapping,
                   by = "locus_tag",
                   all.x = TRUE)

# On rajoute les colonnes du kegg_file pour les pathways 
colnames(kegg_all)<-c("locus_tag","Pathway_ID")
res_annot <- merge(res_annot, kegg_all,
                   by = "locus_tag",
                   all.x = TRUE)
# ===============================



# === Sélection de pathways ===

# Les gènes en lien avec la traduction sont : 
# sao03011 : Ribosome, sao03009 : Ribosome biogenesis, sao03016 : Transfer RNA biogenesis, sao03012 : Translation factors
# On les trouve dans BRITE

translation_genes = !is.na(res_annot$Pathway_ID) &
  (grepl("(^|;)sao03011(;|$)", res_annot$Pathway_ID) |
     grepl("(^|;)sao03009(;|$)", res_annot$Pathway_ID) |
     grepl("(^|;)sao03016(;|$)", res_annot$Pathway_ID) |
     grepl("(^|;)sao03012(;|$)", res_annot$Pathway_ID))


# Filtre uniquement pour les AA_tRNA_synthetases
AA_tRNA_synthetases = !is.na(res_annot$Pathway_ID) &
  grepl("(^|;)sao03016(;|$)", res_annot$Pathway_ID) &
  grepl("-tRNA synthetase", res_annot$product)

res_annot$is_AA_tRNA <- AA_tRNA_synthetases

# Typical genes 
typical_members = subset(
  res_annot,
  symbol %in% c("pth", "infA", "infB", "infC", "frr", "tsf")
)
# ===========================



# ==== VOLCANO PLOT (bonus) ====

# Préparer les colonnes utiles
res_annot$minusLog10Padj <- -log10(res_annot$padj)

res_annot$diffexp <- "Not significant"
res_annot$diffexp[res_annot$log2FoldChange > 0 & res_annot$is_sig] <- "Up"
res_annot$diffexp[res_annot$log2FoldChange < 0 & res_annot$is_sig] <- "Down"

# Top 10 gènes les plus signifs
top10 <- res_annot[order(res_annot$padj), ][1:10, ]

pdf("volcano_repro.pdf",width=8,height=7)
# Plot
volcano <- ggplot(res_annot, aes(x = log2FoldChange, y = minusLog10Padj)) +
  geom_point(aes(color = diffexp), alpha = 0.7, size = 1.5) +
  
  # Seuils
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  # Labels des top gènes
  geom_text_repel(
    data = top10,
    aes(label = symbol),
    size = 2.2,
    segment.color = "black",
    segment.size = 0.8,
    min.segment.length = 0,
    box.padding = 0.8,
    point.padding = 0
  ) +
  
  # Couleurs
  scale_color_manual(
    values = c("Down" = "blue", "Not significant" = "grey70", "Up" = "red")
  ) +
  
  theme_classic() +
  labs(
    title = "Volcano plot : Reproduction results",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"adjusted p-value"),
    color = "Differential expression"
  )+
  theme(plot.title = element_text(hjust = 0.5))

print(volcano)
dev.off()
#======================================



# ===== MA-PLOT TRANSLATION GENES =====
plot_df=res_annot[translation_genes, ]

pdf("MA_plot_translationgenes.pdf", width = 5, height = 5)
p <- ggplot(
  data = plot_df,
  aes(x = log2BaseMean, y = log2FoldChange, color = signif)
) +
  
  # Tous les gènes de traduction (gris / rouge)
  geom_point(alpha = 0.9, size = 1.6) +
  
  # Cercle noir autour des AA-tRNA synthetases (restreint aux translation_genes)
  geom_point(
    data = subset(plot_df, is_AA_tRNA),
    aes(shape = "AA_tRNA_synthetases"),
    color = "black", size = 1.6, stroke = 1.3
  ) +
  
  # Couleurs gris / rouge
  scale_color_manual(values = c("Non-Significant" = "grey60",
                                "Significant" = "red")) +
  
  # Cercle vide pour les AA-tRNA
  scale_shape_manual(values = c("AA_tRNA_synthetases" = 1)) +
  
  # Axes
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-6, 5),
    breaks = seq(-6, 5, by = 1),
    expand = c(0, 0)
  ) +
  
  # Labels des gènes typiques
  geom_text_repel(
    data = typical_members,
    aes(x = log2BaseMean, y = log2FoldChange, label = symbol),
    inherit.aes = FALSE,
    fontface = "italic",
    size = 5,
    segment.color = "black",
    segment.size = 1.2,
    min.segment.length = 0,
    box.padding = 1.6,
    point.padding = 0
  ) +
  
  # Ligne horizontale
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(
    x = expression(log[2]~"BaseMean"),
    y = expression(log[2]~"Fold Change"),
    color = NULL,
    shape = NULL,
    title = ""
  ) +
  
  theme_classic() +
  theme(
    # background supprimé
    panel.background = element_rect(fill = NA, color = NA),
    panel.border     = element_rect(fill = NA, color = "black", size = 1)
  )


# === Gestion de la légende ===

p <- p + theme(
  legend.box = "horizontal"
)

# Ordre des légendes
p <- p + guides(
  color = guide_legend(order = 1),
  shape = guide_legend(order = 2)
)

# Légende en haut (plus lisible)
p <- p + theme(
  legend.position = "top",
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10)
)

print(p)
dev.off()
#======================================

