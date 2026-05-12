process VEP {
    tag "vep_annotation"
    label 'process_vep'
    conda 'bioconda::ensembl-vep=110.1'
    container 'ensemblorg/ensembl-vep:release_110.1'

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
