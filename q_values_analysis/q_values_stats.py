#!/usr/bin/env python3
import os
import gzip
import csv
import urllib.request

# dossier courant = q_values_analysis
WORK_DIR = os.getcwd()
DATA_DIR = os.path.join(WORK_DIR, "fastq_files")
os.makedirs(DATA_DIR, exist_ok=True)

# fichier CSV contenant les URLs
CSV_FILE = os.path.join(WORK_DIR, "sample_url.csv")

with open(CSV_FILE) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        sample_name = row['sample']
        url = row['url']
        gz_fname = os.path.join(DATA_DIR, os.path.basename(url))

        # Télécharger si le fichier n'existe pas déjà
        if not os.path.exists(gz_fname):
            print(f"Téléchargement de {gz_fname} ...")
            urllib.request.urlretrieve(url, gz_fname)

        # Lire le fastq.gz et calculer min, max, mean de Q
        q_scores = []
        with gzip.open(gz_fname, 'rt') as f:
            for i, line in enumerate(f):
                if (i+1) % 4 == 0:  # ligne QUAL
                    q_scores.extend([ord(c) - 33 for c in line.strip()])

        if q_scores:
            q_min = min(q_scores)
            q_max = max(q_scores)
            q_mean = sum(q_scores) / len(q_scores)
            print(f"{sample_name}: min={q_min}, max={q_max}, mean={q_mean:.2f}")
        else:
            print(f"{sample_name}: aucune donnée de qualité trouvée.")
