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


process BOWTIE {
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
    # 1. Indexation du génome (à exécuter une seule fois si l'index n'existe pas)
    # Si l'index n'existe pas :
    # bowtie-build ${ref_fasta} ${index_prefix}
    
    # 2. Alignement avec Bowtie (version 0.12.7)
    # -S = sortie SAM
    bowtie -S -p 4 ${index_prefix} ${fastq_trimmed_file} ${sam_file}

    # 3. Conversion SAM -> BAM non trié (nécessite samtools) car 
    featureCounts prend un BAM en entrée.
    samtools view -bS ${sam_file} > ${bam_unsorted}

    # 4. Tri du fichier BAM (nécessaire pour featureCounts)
    samtools sort ${bam_unsorted} -o ${bam_sorted}
    
    # 5. Nettoyage
    rm ${sam_file} ${bam_unsorted}
    """
    
    // Ajoutez publishDir si vous voulez les BAM dans le répertoire de résultats
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

    container 'mathisrescanieres/featurecount:1.4.6'

    script:
    """
    # Comptage des reads avec featureCounts
    featureCounts \
        -a ${gff} \
        -g ID \
        -o counts_matrix.txt ${bam_files.join(' ')}
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