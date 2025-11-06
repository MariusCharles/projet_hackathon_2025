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

    input:
    file fastq_file

    output:
    path "*.sorted.bam"

    publishDir "results/bam_files", mode: 'copy'

    script:
    """
    set -e
    wget -O genome.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
    gunzip -c genome.fna.gz > genome.fna
    bowtie-build genome.fna genome_index
    bowtie -S genome_index ${fastq_file} > ${fastq_file.simpleName}.sam
    samtools view -bS ${fastq_file.simpleName}.sam > ${fastq_file.simpleName}.bam
    samtools sort ${fastq_file.simpleName}.bam -o ${fastq_file.simpleName}.sorted.bam
    samtools index ${fastq_file.simpleName}.sorted.bam
    rm ${fastq_file.simpleName}.sam ${fastq_file.simpleName}.bam
    """
}


process FEATURECOUNTS {

    input:
    file bam_file                  

    output:
    file "*.txt"

    publishDir "results/coun_matrix", mode: 'copy'
     

    script:
    """
    wget -O genome.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz

    featureCounts -a genome.gff.gz -o ${bam_file.simpleName}.txt ${bam_file}
    """

}

// // Channel pour le fichier GFF
// gff_file = Channel.fromPath('data/GCF_000013425.1_ASM1342v1_genomic.gff')

// process FEATURECOUNTS {

//     input:
//     file bam_files from BOWTIE.out.collect()   
//     file gff from gff_file                      

//     output:
//     file "counts_matrix.txt"

//     container 'mathisrescanieres/featurecount:1.4.6'

//     script:
//     """
//     # Comptage des reads avec featureCounts
//     featureCounts \
//         -a ${gff} \
//         -g ID \
//         -o counts_matrix.txt ${bam_files.join(' ')}
//     """
    
//     publishDir "results", mode: 'copy'
// }


workflow {
sra_url = Channel
    .fromPath('./data/data_URL.txt')
    .flatMap { file -> file.readLines() }  
    .map { it.trim() }                    
fastq_files = DOWNLOAD(sra_url)
fastq_trimmed = CUTADAPT(fastq_files)
bam_files = BOWTIE(fastq_trimmed)
count_txt = FEATURECOUNTS(bam_files)
}
