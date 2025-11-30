# 🧬 Q-values Analysis

Analyse des valeurs de qualité présentes dans des fichiers FASTQ téléchargés depuis SRA. Le but est de déterminer quelles valeurs de Q sont réellement utilisées dans le pipeline RNA-seq de l'article.

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

Le dossier `fastq_files/` nest pas présent par défaut dans le dossier `q_values_analysis` mais est créé une fois les scripts Python lancés.

## Données d’entrée : sample_url.csv

Le fichier CSV contient deux colonnes : `sample` et `url`.  
Les scripts utilisent ces URL pour télécharger automatiquement les FASTQ.gz.  

## Scripts disponibles

### `q_unique_values.py` : Valeurs uniques des Q-scores

Ce script :

- télécharge les FASTQ.gz (via `wget`)

- parcourt uniquement les lignes de qualité (1 ligne sur 4)

- convertit les caractères ASCII en Phred scores (ord(c) - 33)

- collecte les valeurs Q uniques par fichier d'entrée

**Résultats (q_unique_values_results.txt)**

Tous les fichiers `.fastq` contiennent les mêmes valeurs uniques de Q : 6 valeurs Q uniques → [2, 14, 22, 27, 33, 37]

### `q_values_stats.py` → Statistiques globales (min / max / mean)

Ce script : 

- télécharge les `fastq.gz` (via `wget`)

- parcourt uniquement les lignes de qualité (1 ligne sur 4)

- convertit les caractères ASCII en Phred scores (ord(c) - 33)

- calcule les statistiques des valeurs de Q par fichier d'entrée

**Interprétation**:

- max global = 37 → très haute qualité, proche du plafond Illumina.

- min global = 2 → très rare : probablement quelques bases très faibles.

- mean globale ≈ 35.5 → qualité globalement excellente.

## Valeurs de Q-score à intégrer au pipeline

```
[1, 3, 15, 23, 28, 34]
```

Ces valeurs couvrent toutes les zones critiques de la distribution réelle des Q-scores observés dans les FASTQ.  
Chaque valeur unique testée est plus ou moins permissive. Ainsi :

| Q testé |  Intérêt |
|--------|---------|
| **1**  | Vérifie l’impact d’un trimming ultra permissif |
| **3**  | Test raisonnable bas niveau |
| **15** | Cas médian, trimming modéré |
| **23** | Trimming plus strict |
| **28** | Très strict |
| **34** | Ultra strict (risque d’élaguage élevé) |


## 🔁 Reproductibilité totale

Dans un souci de reproductibilité complète de cette analyse, il est recommander de créer un environnement `conda` dédié à l'exécution des scripts python.

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

### Rendre les scripts exécutables

Après avoir cloné le dépôt, Se placer dans le dossier `q_values_analysis/` contenant les scripts. Appliquer les permissions d’exécution à tous les scripts Python :

```bash
chmod +x *.py
```

### Lancer les scripts

Toujours dans le dossier `q_values_analysis/` des les scripts :

```bash
./q_unique_values.py
./q_values_stats.py
```

Les FASTQ seront téléchargés automatiquement dans le dossier `fastq_files/` et les résultats imprimés à l’écran.
