# ===========================================================
# Dockerfile: featureCounts 1.4.6-p3 (Subread) sur Ubuntu 14.04
# ===========================================================

FROM ubuntu:14.04

# Informations de l'image
LABEL creator="mathis.rescanieres"
LABEL description="Image Docker contenant featureCounts 1.4.6-p3 (Subread package)"

# Eviter les blocages lors de l'installation
ENV DEBIAN_FRONTEND=noninteractive

# Installer les dépendances nécessaires
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    tar \
    && rm -rf /var/lib/apt/lists/*

# Se placer dans /opt pour télécharger et compiler
WORKDIR /opt

# Télécharger, compiler et installer featureCounts
RUN wget -O subread.tar.gz "https://downloads.sourceforge.net/project/subread/subread-1.4.6-p3/subread-1.4.6-p3-source.tar.gz" && \
    tar -xzf subread.tar.gz && \
    cd subread-1.4.6-p3-source/src && \
    make -f Makefile.Linux featureCounts && \
    mv featureCounts /usr/local/bin/ && \
    cd /opt && rm -rf subread-1.4.6-p3-source subread.tar.gz

# Vérifier l'installation
RUN featureCounts -v

# Répertoire de travail pour les données
WORKDIR /data

# Point d'entrée par défaut
ENTRYPOINT ["featureCounts"]
