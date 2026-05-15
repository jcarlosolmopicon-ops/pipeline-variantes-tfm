process VEP {
    tag "vep_annotation"
    label 'process_vep'
    publishDir "${params.outdir}/vep", mode: 'copy'

    conda 'bioconda::ensembl-vep=115'

    input:
    path vcf

    output:
    path "annotated.vcf.gz",              emit: vcf
    path "annotated.vcf.gz_summary.html", emit: summary

    script:
    """
    vep \\
        --input_file ${vcf} \\
        --output_file annotated.vcf.gz \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --assembly GRCh38 \\
        --cache \\
        --dir_cache \$HOME/.vep \\
        --species homo_sapiens \\
        --everything \\
        --fork ${task.cpus}
    """
}
