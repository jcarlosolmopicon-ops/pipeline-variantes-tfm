process BQSR {
    tag "${meta.id}"
    label 'process_gatk'
    publishDir "${params.outdir}/bqsr", mode: 'copy'

    conda 'bioconda::gatk4=4.5.0.0 bioconda::samtools=1.19'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path genome
    path genome_index
    path dbsnp
    path dbsnp_index

    output:
    tuple val(meta), path("${meta.id}.bqsr.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.bqsr.bam.bai"), emit: bai
    path "${meta.id}.recal.table",                     emit: table

    script:
    """
    gatk BaseRecalibrator \\
        -I ${bam} \\
        -R ${genome} \\
        --known-sites ${dbsnp} \\
        -O ${meta.id}.recal.table

    gatk ApplyBQSR \\
        -I ${bam} \\
        -R ${genome} \\
        --bqsr-recal-file ${meta.id}.recal.table \\
        -O ${meta.id}.bqsr.bam

    samtools index ${meta.id}.bqsr.bam
    """
}
