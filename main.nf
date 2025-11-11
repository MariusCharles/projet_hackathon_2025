#!/usr/bin/env nextflow

process DOWNLOAD_FASTQ {
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

process DOWNLOAD_GENOME {
    output:
    file "genome.fna"

    script:
    """
    set -e
    wget -O genome.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
    gunzip -c genome.fna.gz > genome.fna
    """
}

process BUILD_INDEX {
    input:
    file genome_fna

    output:
    path "genome_index.*"

    script:
    """
    bowtie-build ${genome_fna} genome_index
    """
}


process BOWTIE {

    input:
    tuple file(fastq_file), file(genome_index_files)

    output:
    path "*.sorted.bam"
    path "*.sorted.bam.bai"

    publishDir "results/bam_files", mode: 'copy'

    script:
    """
    set -e
    bowtie -S genome_index ${fastq_file} > ${fastq_file.simpleName}.sam
    samtools view -bS ${fastq_file.simpleName}.sam > ${fastq_file.simpleName}.bam
    samtools sort ${fastq_file.simpleName}.bam -o ${fastq_file.simpleName}.sorted.bam
    samtools index ${fastq_file.simpleName}.sorted.bam
    rm ${fastq_file.simpleName}.sam ${fastq_file.simpleName}.bam
    """
}


process DOWNLOAD_GTF {

    output:
    path "genome.gtf"

    script:
    """
    set -e
    wget -O genome.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gtf.gz
    gunzip genome.gtf.gz
    """
}


process FEATURECOUNTS {

    input:
    tuple val(bam_list), path(gtf_file)                

    output:
    file "counts_matrix.txt"

    publishDir "results/count_matrix", mode: 'copy'
     
    script:
    """
    set -e
    BAMS=$(echo ${bam_list.join(' ')})
    featureCounts -a ${gtf_file} -t gene -g ID -s 1 -o counts_matrix.txt $BAMS
    """
}

process DESEQ {
    publishDir 'results/deseq2', mode: 'copy'

    input:
    file counts_matrix
    file coldata_csv
    file deseq_script

    output:
    file 'deseq2_results.csv'
    file 'MA_plot.pdf'

    script:
    """
      eval "\$(micromamba shell hook -s bash)"
      micromamba activate base

      Rscript ${deseq_script} ${counts_matrix} ${coldata_csv}
    """
}




workflow {
// Get les urls dans config.csv
sra_url = Channel
    .fromPath('./data/config.csv')
    .splitCsv(header: true)
    .map { row -> row.url.trim() }                
// On télécharge les fastq (6 channels)
fastq_files = DOWNLOAD_FASTQ(sra_url)
// On trim les fastq (6 channels)
fastq_trimmed = CUTADAPT(fastq_files)
// On télécharge le génome (fasta) 1 seule fois
genome_fna = DOWNLOAD_GENOME()
// On build l'index du génome 1 seule fois
genome_index = BUILD_INDEX(genome_fna)
// On associe les fastq aux genome_index 
reads_with_index = fastq_trimmed.combine(genome_index)
// On aligne les fastq => bam 
bam_files = BOWTIE(reads_with_index)
// On fait collect => on attend l'output de toutes les channels et all_bams contient tous les bams 
all_bams = bam_files.collect()
// On download le GTF
genome_gtf = DOWNLOAD_GTF()
// On associe bams + gtf
bam_with_gtf = all_bams.combine(genome_gtf)
// On génère la matrice de comptes 
count_txt = FEATURECOUNTS(bam_with_gtf)
// On lit la matrice de comptes + coldata et on run Deseq + plot 
deseq_results = DESEQ(
    count_txt,                 
    file('./data/config.csv'), 
    file('run_deseq.R')        // pour l'instant script dans le répertoire courant
)
}

