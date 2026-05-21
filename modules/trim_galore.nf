process TRIM_GALORE {
    tag "${meta.id}"
    label 'low_cpu'
    publishDir "${params.outdir}/trim", mode: 'copy'
    conda 'bioconda::trim-galore=0.6.10'
    errorStrategy 'ignore'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("*_val_{1,2}.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt",                emit: trim_log
    script:
    """
    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --gzip \\
        ${reads}
    """
}
