process VEP {
    tag "vep_annotation"
    label 'vep'
    publishDir "${params.outdir}/vep", mode: 'copy'
    conda 'bioconda::ensembl-vep=110.1'
    container 'ensemblorg/ensembl-vep:release_110.1'
    errorStrategy 'ignore'
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
        --offline \\
        --dir_cache /home/olmop/.vep \\
        --species homo_sapiens \\
        --sift b \\
        --polyphen b \\
        --symbol \\
        --numbers \\
        --biotype \\
        --canonical \\
        --variant_class \\
        --fork ${task.cpus}
    """
}
