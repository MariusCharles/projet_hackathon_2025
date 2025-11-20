FROM ubuntu:18.04

# Installer les dépendances de base
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    samtools \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Télécharger et installer les binaires précompilés de Bowtie 0.12.7
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download -O bow>
    unzip bowtie-0.12.7-linux-x86_64.zip && \
    cp bowtie-0.12.7/bowtie* /usr/local/bin/ && \
    chmod +x /usr/local/bin/bowtie*

# Vérification
RUN bowtie --version