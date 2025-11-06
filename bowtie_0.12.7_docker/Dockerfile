# 1. Utilise une image de base Ubuntu ancienne pour la compatibilité
FROM ubuntu:14.04

# 2. Installe les dépendances nécessaires à la compilation (make, gcc, zlib)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    make \
    #gcc \
    #g++ \
    unzip \
    zlib1g-dev gcc-4.4 g++-4.4\
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Définit Bash comme shell par défaut pour toutes les commandes RUN suivantes
SHELL ["/bin/bash", "-c"]

# 4. Télécharge et décompresse la version spécifique de Bowtie
WORKDIR /opt/
RUN wget http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download -O bowtie-0.12.7-src.zip
RUN unzip bowtie-0.12.7-src.zip
RUN rm bowtie-0.12.7-src.zip


# 5. Télécharge et installe Samtools nécessaire pour convertir sortie SAM de bowtie en BAM
#utilisable par FeatureCount
WORKDIR /opt/
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
RUN tar -xjf samtools-1.3.1.tar.bz2
WORKDIR /opt/samtools-1.3.1
# Samtools a besoin des dépendances de compilation et de libncurses 
RUN make && make install
RUN rm /opt/samtools-1.3.1.tar.bz2

# 6. Compile Bowtie à partir des sources
WORKDIR /opt/bowtie-0.12.7
# gérer un bug connu de bowtie avec les versions récentes de GCC.
RUN make CC=gcc-4.4 CXX=g++-4.4   

# 7. Définit le chemin d'accès à l'exécutable
ENV PATH="/opt/bowtie-0.12.7:$PATH"

# 8. Commande par défaut (optionnel)
CMD ["bash", "-c", "bowtie", "--version"]
