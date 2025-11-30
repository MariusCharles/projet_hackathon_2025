#!/usr/bin/env python3
import os
import gzip
import csv
import urllib.request
import subprocess

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
            subprocess.run(["wget", "-c", url, "-O", gz_fname], check=True)

        # Lire le fastq.gz et collecter les valeurs Q uniques
        unique_q = set()
        with gzip.open(gz_fname, 'rt') as f:
            for i, line in enumerate(f):
                if (i + 1) % 4 == 0:  # ligne QUAL
                    unique_q.update([ord(c) - 33 for c in line.strip()])

        if unique_q:
            sorted_q = sorted(unique_q)
            print(f"{sample_name}: {len(unique_q)} valeurs Q uniques → {sorted_q}")
        else:
            print(f"{sample_name}: aucune donnée de qualité trouvée.")
