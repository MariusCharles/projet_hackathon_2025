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

    cpus { Math.max(1, (params.core_nb / 4).toInteger()) }

    input:
    file fastq_file 
    val q  

    output:
    file "*.fastq"

    publishDir "results/fastq_files_trimmed", mode: 'copy'

    script:
    """
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
             -m 25 -q ${q} \
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

    cpus { Math.max(1, (params.core_nb / 2).toInteger()) }

    input:
    path fastq_file
    path genome_index

    output:
    path "*.sorted.bam", emit: sorted_bam
    path "*.bai", emit: bai_index


    publishDir "results/bam_files", mode: 'copy'

    script:
    """
    set -e
    bowtie -S genome_index ${fastq_file} | \
    samtools view -bS - | \
    samtools sort -@ ${task.cpus} -o ${fastq_file.simpleName}.sorted.bam -
    samtools index ${fastq_file.simpleName}.sorted.bam
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
    
    cpus { Math.max(1, (params.core_nb / 2).toInteger()) }

    input:
    path bam_files
    path gtf_file

    output:
    file "counts_matrix.txt"

    publishDir "results/count_matrix", mode: 'copy'

    script:
    """
    featureCounts -T ${task.cpus} \
        -F GTF -a ${gtf_file} -t gene -g gene_id -s 1 \
        -o counts_matrix.txt ${bam_files.join(' ')}
    """
}


process DESEQ {
    publishDir 'results/deseq2', mode: 'copy'

    input:
    file counts_matrix
    file coldata_csv

    output:
    file 'deseq2_results.csv'
    file 'MA_plot.pdf'

    script:
    """
    eval "\$(micromamba shell hook -s bash)"
    micromamba activate base

    bin/run_deseq.R ${counts_matrix} ${coldata_csv}
    """
}



workflow {
// Get les urls dans config.csv
sra_url = Channel
    .fromPath(params.config)
    .splitCsv(header: true)
    .map { row -> row.url.trim() }      

// On télécharge les fastq (6 channels)
fastq_files = DOWNLOAD_FASTQ(sra_url)

// On trim les fastq (6 channels)
fastq_trimmed = CUTADAPT(fastq_files, params.q)

// On télécharge le génome (fasta) 1 seule fois
genome_fna = DOWNLOAD_GENOME()

// On build l'index du génome 1 seule fois
genome_index = BUILD_INDEX(genome_fna)

// On aligne les fastq => bam 
sorted_bams = BOWTIE(fastq_trimmed,genome_index).sorted_bam

// On fait collect => on attend que tous les BAM soient produits
all_bams = sorted_bams.collect()

// On download le GTF
genome_gtf = DOWNLOAD_GTF()

// On génère la matrice de comptes 
count_txt = FEATURECOUNTS(all_bams,genome_gtf)

// On lit la matrice de comptes + coldata et on run Deseq + plot 
deseq_results = DESEQ(count_txt, file(params.config))
}
