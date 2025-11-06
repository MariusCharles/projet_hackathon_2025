# Utiliser une image légère avec Debian
FROM debian:bullseye-slim

# Installer wget
RUN apt-get update && \
    apt-get install -y wget ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Définir le dossier de travail par défaut
WORKDIR /data

# Définir le point d'entrée comme wget
ENTRYPOINT ["wget"]

