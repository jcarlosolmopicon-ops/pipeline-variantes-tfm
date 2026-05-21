process FASTQC {
    tag "${meta.id}"
    label 'high_cpu'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    conda 'bioconda::fastqc=0.12.1'
    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

script:
"""
export NXF_DEBUG=\${NXF_DEBUG:-0}
fastqc \\
    --threads ${task.cpus} \\
    --outdir . \\
    ${reads}
"""
}
