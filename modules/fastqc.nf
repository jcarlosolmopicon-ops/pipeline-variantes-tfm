process FASTQC {
    tag "${meta.id}"
    label 'high_cpu'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads}
    """
}
