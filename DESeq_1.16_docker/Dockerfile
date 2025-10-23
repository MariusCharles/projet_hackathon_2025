FROM mambaorg/micromamba:1.5.8 

# Installer R 3.4.1 + DESeq2 1.16 + ggplot2
# + runtimes (libgfortran.so.3, libgcc, libstdc++)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    r-base=3.4.1 \
    bioconductor-deseq2=1.16.* \
    r-ggplot2=3.0.* \
    libgfortran=3.0.0 \
    libgcc-ng=7.* \
    libstdcxx-ng=7.* \
 && micromamba clean -a -y

# Pour que les commandes suivantes utilisent bash + conda env "base"
SHELL ["/bin/bash", "-lc"]
