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

workflow {
sra_ids = Channel
    .fromPath('data_SRA.txt')
    .flatMap { file -> file.readLines() }  
    .map { it.trim() }                    
fastq_files = FASTQ_DOWN(sra_ids)
fastq_trimmed = CUTADAPT(fastq_files)
}