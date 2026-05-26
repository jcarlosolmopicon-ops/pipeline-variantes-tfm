process MARKDUPLICATES {
    tag "${meta.id}"
    label 'gatk'
    publishDir "${params.outdir}/markdup", mode: 'symlink'
    conda 'bioconda::gatk4=4.5.0.0 bioconda::samtools=1.19'
    container 'broadinstitute/gatk:4.5.0.0'
    errorStrategy 'ignore'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path("${meta.id}.markdup.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.markdup.bam.bai"), emit: bai
    path "${meta.id}.markdup.metrics",                    emit: metrics
    script:
    """
    gatk --java-options "-Xmx16g" MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.markdup.bam \\
        -M ${meta.id}.markdup.metrics \\
        --TMP_DIR .
    samtools index ${meta.id}.markdup.bam
    """
}
