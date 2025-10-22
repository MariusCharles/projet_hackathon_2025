#! /usr/bin/env nextflow

params.input = "data/*_R{1,2}.fastq"
params.outdir = "results"

channel
    .fromFilePairs(params.input,flat=true)
    .set{paired_reads}

process cutadapt{
    tag {sample_id}

    container 'quay.io/biocontainers/cutadapt:1.11-??'

    input:
    tuple val(sample_id),path(reads)

    output:
    tuple val(sample_id),path(["nettoye_R1.fastq", "nettoye_R2.fastq"])

    script:
    """
    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o nettoye_R1.fastq -p nettoye_R2.fastq \
    -m 25 -q ?? \
    ${reads[0]} ${reads[1]}
    """

    publishDir "results/nettoye_etape_1", mode: 'copy'
}