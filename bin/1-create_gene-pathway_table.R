library(KEGGREST) 

BASE_PATH = "/Users/mafaldafrere/Documents/Cours/IODAA/HACKATHON/PROJET/clean_KEGG"

# On récupère les associations gènes - pathways 
gene_pathway= keggLink(target = "brite", source = "sao") 

# On retire le préfixe 'sao'
Gene_ID = sub(pattern = "^sao:",  replacement = "", names(gene_pathway)) 
# On retire le préfixe 'br'
Pathway_ID = sub(pattern = "^br:",  replacement = "", gene_pathway) 

# To dataframe
gene_pathway= data.frame(Gene_ID, Pathway_ID) 

# On garde un seul gène par ligne et on lui associe ses différents pathways 
gene_pathway = aggregate(. ~ Gene_ID, data = gene_pathway,
                            FUN = function(x) paste(unique(x), collapse = ";"))

# Save la table obtenue
write.table(gene_pathway, file = file.path(BASE_PATH,"gene_pathway_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

