#!/usr/bin/env nextflow

process DOWNLOAD {
    input:
    val url

    output:
    file "*.fastq.gz"  

    publishDir "results/fichier_fastq", mode: 'copy'

    script:
    """
    wget -c ${url} -P .
    """
}


process CUTADAPT {

    input:
    file fastq_file   

    output:
    file "*.fastq"

    publishDir "results/fastq_files_trimmed", mode: 'copy'

    script:
    """
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
             -m 25 -q 20 \
             -o ${fastq_file.simpleName}_trimmed.fastq ${fastq_file}
    """
    
}


process BOWTIE {

    tag "${fastq_trimmed_file.baseName}"

    input:
    file fastq_trimmed_file  //fichier FASTQ trimé
    file ref_fasta    //fichier FASTA de référence
    val index_prefix   //prefixe de l'index Bowtie

    output:
    file "${fastq_trimmed_file.baseName}.bam"

    script:
    def sam_file = "${fastq_trimmed_file.baseName}.sam"
    def bam_unsorted = "${fastq_trimmed_file.baseName}_unsorted.bam"
    def bam_sorted = "${fastq_trimmed_file.baseName}.bam"

    """
    # 1️⃣ Création de l'index Bowtie si nécessaire
    if [ ! -f "${index_prefix}.1.ebwt" ]; then
        bowtie-build ${ref_fasta} ${index_prefix}
    fi

    # 2️⃣ Alignement FASTQ → SAM
    bowtie -S -p 4 ${index_prefix} ${fastq_trimmed_file} > ${sam_file}

    # 3️⃣ Conversion SAM → BAM non trié
    samtools view -bS ${sam_file} -o ${bam_unsorted}

    # 4️⃣ Tri du BAM
    samtools sort ${bam_unsorted} -o ${bam_sorted}

    # 5️⃣ Indexation du BAM trié
    samtools index ${bam_sorted}

    # 6️⃣ Vérification du contenu du BAM (affiche les premières lignes)
    samtools view ${bam_sorted} | head

    # 7️⃣ Affichage des statistiques de base
    samtools flagstat ${bam_sorted}

    # 8️⃣ Nettoyage des fichiers intermédiaires
    rm ${sam_file} ${bam_unsorted}
    """

    publishDir "results", mode: 'copy'
}

// Channel pour le fichier GFF
gff_file = Channel.fromPath('data/GCF_000013425.1_ASM1342v1_genomic.gff')

process FEATURECOUNTS {

    input:
    file bam_files from BOWTIE.out.collect()   
    file gff from gff_file                      

    output:
    file "counts_matrix.txt"

    publishDir "results", mode: 'copy'

    script:
    """
    # Comptage des reads avec featureCounts
    featureCounts \
        -a ${gff} \
        -g ID \
        -o counts_matrix.txt ${bam_files.join(' ')}
    """
}



workflow {
sra_url = Channel
    .fromPath('./data/data_url.txt')
    .flatMap { file -> file.readLines() }  
    .map { it.trim() }                    
fastq_files = DOWNLOAD(sra_url)
fastq_trimmed = CUTADAPT(fastq_files)
// bam_files = ALIGNMENT(fastq_trimmed)   // Partie de Eliott
counts = FEATURECOUNTS(bam_files)
}