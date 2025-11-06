#!/usr/bin/env nextflow

process FASTQ_DOWN {
    input:
    val sra_id

    output:
    file "${sra_id}_*.fastq"  

    script:
    """
    fasterq-dump ${sra_id} -O .
    """
}


process CUTADAPT {

    input:
    file fastq_file   

    output:
    file "${fastq_file.baseName}_trimmed.fastq"

    script:
    """
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
             -m 25 -q 20 \
             -o ${fastq_file.baseName}_trimmed.fastq ${fastq_file}
    """

    publishDir "results", mode: 'copy'
}

process FEATURECOUNTS {

    input:
    file bam_file

    output:
    file "${bam_file.baseName}_counts.txt"

    script:
    """
    featureCounts -a annotation.gtf -o ${bam_file.baseName}_counts.txt ${bam_files.join(' ')}
    """
    
    publishDir "results", mode: 'copy'
}


workflow {
sra_ids = Channel
    .fromPath('data_SRA.txt')
    .flatMap { file -> file.readLines() }  
    .map { it.trim() }                    
fastq_files = FASTQ_DOWN(sra_ids)
fastq_trimmed = CUTADAPT(fastq_files)
// bam_files = ALIGNMENT(fastq_trimmed)   // Partie de Eliott
counts = FEATURECOUNTS(bam_files)
}