# Dockerfile pour Bowtie version 0.12.7 et Samtools 1.3.1
# ----------------------------------------------------

# 1. Utilise une image de base Ubuntu ancienne pour la compatibilité
FROM ubuntu:14.04

# 2. Installe les dépendances nécessaires à la compilation (make, gcc, zlib)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    make \
    unzip \
    zlib1g-dev \
    # ASSURE LA PRÉSENCE DE LA COMMANDE GÉNÉRIQUE 'gcc' et 'g++'
    gcc g++ \
    # Garde les versions 4.4 pour la compilation spécifique de Bowtie
    gcc-4.4 g++-4.4 \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Définit Bash comme shell par défaut pour toutes les commandes RUN suivantes
SHELL ["/bin/bash", "-c"]

# --- AJOUT DE SAMTOOLS ---
# 4. Télécharge et installe Samtools
WORKDIR /opt/
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
RUN tar -xjf samtools-1.3.1.tar.bz2
WORKDIR /opt/samtools-1.3.1
# La commande 'make' utilise maintenant le 'gcc' générique installé à l'étape 2.
RUN make && make install
RUN rm /opt/samtools-1.3.1.tar.bz2
# --- FIN DE L'AJOUT SAMTOOLS ---

# 5. Télécharge et décompresse la version spécifique de Bowtie
WORKDIR /opt/
RUN wget http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download -O bowtie-0.12.7-src.zip
RUN unzip bowtie-0.12.7-src.zip
RUN rm bowtie-0.12.7-src.zip

# 6. Compile Bowtie à partir des sources
WORKDIR /opt/bowtie-0.12.7
# Utilise explicitement les versions 4.4 pour gérer le bug connu de Bowtie.
RUN make CC=gcc-4.4 CXX=g++-4.4   

# 7. Définit le chemin d'accès aux exécutables (Bowtie et Samtools)
ENV PATH="/opt/bowtie-0.12.7:/usr/local/bin:$PATH"

# 8. Commande par défaut (optionnel)
CMD ["bash", "-c", "bowtie --version && samtools --version"]
