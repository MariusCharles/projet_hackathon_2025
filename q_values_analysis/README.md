# 🧬 Q-values Analysis

Analyse des valeurs de qualité Phred présentes dans des fichiers FASTQ téléchargés depuis SRA afin de déterminer quelles valeurs de Q sont réellement utilisées dans le pipeline RNA-seq de l'article.

## Contenu du dépôt

```bash
.
├── q_unique_values.py            # Extraction des valeurs Q uniques par FASTQ
├── q_unique_values_results.txt   # Résultats associés
├── q_values_stats.py             # min / max / mean des Phred scores
├── q_values_stats_results.txt    # Résultats associés
├── sample_url.csv                # Liste des samples et URLs (ENA/SRA)
├── fastq_files/                  # FASTQ.gz téléchargés automatiquement
└── README.md                     # Ce fichier
```

## Données d’entrée : sample_url.csv

Le fichier CSV contient deux colonnes : `sample` et `url`.  
Les scripts utilisent ces URL pour télécharger automatiquement les FASTQ.gz.  

## Scripts disponibles

### q_unique_values.py → Valeurs Q uniques

Ce script :

- télécharge les FASTQ.gz (via `wget`)

- parcourt uniquement les lignes de qualité (1 ligne sur 4)

- convertit les caractères ASCII en Phred scores (ord(c) - 33)

- collecte les valeurs Q uniques par fichier d'entrée

**Résultats (q_unique_values_results.txt)**

Tous les fichiers `.fastq` contiennent les mêmes valeurs uniques de Q : 6 valeurs Q uniques → [2, 14, 22, 27, 33, 37]

### q_values_stats.py → Statistiques globales (min / max / mean)

Ce script : 

- télécharge les FASTQ.gz (via `wget`)

- parcourt uniquement les lignes de qualité (1 ligne sur 4)

- convertit les caractères ASCII en Phred scores (ord(c) - 33)

- calcule les statistiques des valeurs de Q par fichier d'entrée

**Interprétation rapide**:

- max = 37 → très haute qualité, proche du plafond Illumina.

- min = 2 → très rare : probablement quelques bases très faibles.

- mean ≈ 35.6 → qualité globalement excellente.

### Valeurs de Q à intégrer au pipeline

```
[1, 3, 15, 23, 28, 34]
```

Ces valeurs couvrent toutes les zones critiques de la distribution réelle des Q-scores observés dans les FASTQ.  
Elles correspondent aux intervalles suivants :

| Q testé | Correspondance avec les Q trouvés | Intérêt |
|--------|-----------------------------------|---------|
| **1**  | < min réel                        | Vérifie l’impact d’un trimming ultra permissif |
| **3**  | proche du min réel (=2)           | Test raisonnable bas niveau |
| **15** | entre 14 et 22                    | Cas médian, trimming modéré |
| **23** | entre 22 et 27                    | Trimming plus strict |
| **28** | entre 27 et 33                    | Très strict |
| **34** | entre 33 et 37                    | Ultra strict (risque d’élaguage élevé) |


## 🔁 Reproductibilité totale

### Création d'un environnement Conda

```bash
conda create -n qvalues_env python=3.12 -y
conda activate qvalues_env
```

### Téléchargement des dépendances nécessaires

Tout est dans la bibliothèque standard sauf `wget`.

```bash
sudo apt-get update
sudo apt-get install wget
```

### Lancer les scripts

Se placer dans le dossier `q_values_analysis/` contenant les scripts :

```bash
python q_unique_values.py
```

puis 

```bash
python q_values_stats.py
```

Les FASTQ seront téléchargés automatiquement dans le dossier `fastq_files/` et les résultats imprimés à l’écran.