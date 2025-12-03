# 🧬 Projet Hackathon 2025 

Ce projet reproduit le pipeline **RNA-seq** de l’article *“Intracellular Staphylococcus aureus persisters upon antibiotic exposure”* de manière **reproductible**, depuis les fichiers FASTQ jusqu’aux analyses différentielles (MA-plots DESeq2).

L’objectif est de garantir la reproductibilité complète des résultats : mêmes versions de packages, mêmes outils.

---

## En cours

- [ ] Finir le rapport (ajouter résultats, commenter)
- [ ] Commencer le diapo et présentation orale
- [ ] Faire tourner le workflow pour plusieurs valeurs de q

---

## Structure du dépôt

| Dossier / fichier | Description |
|-------------------|--------------|
| **`main.nf`** | Script **Nextflow** principal décrivant l’ensemble du pipeline. |
| **`nextflow.config`** | Fichier de configuration du pipeline : définit les images Docker à utiliser (appel depuis DockerHub). |
| **`docker/`** | Contient un Dockerfile pour chaque outil utilisé (Bowtie, Cutadapt, FeatureCounts, DESeq2, etc.). Les images correspondantes sont disponibles sur DockerHub. |
| **`data/`** | Données d’entrée :<br>• `config.csv` - table de description des échantillons (nom, URL FASTQ, réplicat, condition). |
| **`bin/`** | Scripts exécutables du pipeline :<br>• `1-create_gene-pathway_table.R` : récupération et formatage des associations gènes–pathways KEGG<br>• `2-deseq_table_suppplot.R` : analyse différentielle avec DESeq2 + génération du MA-plot global<br>• `3-create_downstream_plots.R` : génération des volcano plots et MA-plots ciblés<br>• `4-paper_results_comp.R` : comparaison avec les résultats du papier (Venn diagrams, corrélations, scatter plots). |
| **`results_all_q_values/`** | Résultats des runs du pipeline (plots et tables) pour les différentes valeurs de Q-score possibles. |

---

## Objectif

Reproduire les figures principales de l’article à partir des données publiques, dans un environnement **conteneurisé** et **traçable** via Nextflow + Docker.

---

## Exécution rapide (exemple)

### Paramètres de la machine virtuelle utilisée
- **16 CPU**  
- **64 Go de RAM**  
- **400 Go de stockage**

### Installer Git sur la VM (si nécessaire)

Si Git n'est pas installé sur la VM, utiliser la commande suivante pour l’installer :

```bash
sudo apt update
sudo apt install git -y
```

### Cloner le repository Git dans la VM
```bash
git clone https://github.com/MariusCharles/projet_hackathon_2025.git
```

### Rentrer dans le dossier du projet
```bash
cd projet_hackathon_2025
```

### Activer Nextflow dans la VM avec Conda 
```bash
conda init
source ~/.bashrc
conda activate nextflow
```

### Lancer le pipeline

```bash
nextflow run main.nf
```
