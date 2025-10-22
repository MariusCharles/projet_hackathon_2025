# Dockerfile pour Bowtie version 0.12.7
# ------------------------------------

# 1. Utilise une image de base Ubuntu ancienne pour la compatibilité
FROM ubuntu:14.04

# 2. Installe les dépendances nécessaires à la compilation (make, gcc, zlib)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    make \
    gcc \
    g++ \
    unzip \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Télécharge et décompresse la version spécifique de Bowtie
WORKDIR /opt/
RUN wget http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download -O bowtie-0.12.7-src.zip
RUN unzip bowtie-0.12.7-src.zip
RUN rm bowtie-0.12.7-src.zip

# 4. Compile Bowtie à partir des sources
WORKDIR /opt/bowtie-0.12.7
RUN make

# 5. Définit le chemin d'accès à l'exécutable
ENV PATH="/opt/bowtie-0.12.7:$PATH"

# 6. Commande par défaut (optionnel)
CMD ["bowtie"]
