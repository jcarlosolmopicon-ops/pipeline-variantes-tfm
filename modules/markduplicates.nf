process MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_gatk'
    conda 'bioconda::gatk4=4.5.0.0 bioconda::samtools=1.19'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.markdup.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.markdup.bam.bai"), emit: bai
    path "${meta.id}.markdup.metrics",                    emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.markdup.bam \\
        -M ${meta.id}.markdup.metrics \\
        --TMP_DIR /tmp

    samtools index ${meta.id}.markdup.bam
    """
}
