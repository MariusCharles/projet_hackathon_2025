# Dockerfile pour Cutadapt 1.11
FROM python:3.6-slim

# Installer les dépendances nécessaires
RUN apt-get update && \
    apt-get install -y build-essential zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Installer Cutadapt 1.11 via pip
RUN pip install cutadapt==1.11

# Ajouter Cutadapt au PATH (déjà fait par pip)
ENV PATH="/usr/local/bin:$PATH"

# Définir le point d'entrée
ENTRYPOINT ["cutadapt"]
