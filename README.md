# 🧬 Projet Hackathon 2025 

Ce projet reproduit le pipeline **RNA-seq** de l’article *“Intracellular Staphylococcus aureus persisters upon antibiotic exposure”* de manière **reproductible**, depuis les fichiers FASTQ jusqu’aux analyses différentielles (MA-plots DESeq2).

L’objectif est de garantir la reproductibilité complète des résultats : mêmes versions de packages, mêmes outils.

---

## En cours

- [ ] Faire le diapo pour séance 14/11
- [ ] Automatiser nb de cpu utilisé par chaque process (task.cpus...)
- [ ] Construire base de données KEGG translation pour le plot annoté 
- [ ] Gérer conversions noms de gènes
- [ ] Nettoyer le dossier `to_delete/` après validation  
- [ ] Commencer à écrire le rapport

---

## Structure du dépôt

| Dossier / fichier | Description |
|-------------------|--------------|
| **`main.nf`** | Script **Nextflow** principal décrivant l’ensemble du pipeline. |
| **`nextflow.config`** | Fichier de configuration du pipeline : définit les images Docker à utiliser (appel depuis DockerHub). |
| **`docker/`** | Contient un Dockerfile pour chaque outil utilisé (Bowtie, Cutadapt, FeatureCounts, DESeq2, etc.). Les images correspondantes sont disponibles sur DockerHub. |
| **`data/`** | Données d’entrée :<br>• `config.csv` - table de description des échantillons (nom, URL FASTQ, réplicat, condition)<br>• script R (analyse DESeq2). |
| **`results/`** | Résultats du dernier run du pipeline : matrice de comptages, résultats DESeq2, MA-plot. |
| **`to_delete/`** | Dossier temporaire pour fichiers/données à valider avant suppression définitive. |

---

## Objectif

Reproduire les figures principales de l’article à partir des données publiques, dans un environnement **conteneurisé** et **traçable** via Nextflow + Docker.

---

## Exécution rapide (exemple)

### Paramètres de la machine virtuelle utilisée
- **16 CPU**  
- **64 Go de RAM**  
- **400 Go de stockage**

### Fichiers nécessaires
- **`main.nf`**, **`nextflow.config`**, **`data/`**

### Lancer le pipeline

```bash
nextflow run main.nf
