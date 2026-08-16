process CLASSIFY {
    tag "ml_classification"
    label 'classify'
    publishDir "${params.outdir}/classify", mode: 'copy'
    conda "${projectDir}/entorno/environment_ubuntu_2026-05.yml"
    input:
    path vcf
    path model_dir
    output:
    path "hcc1395_annotated_with_ml.tsv", emit: tsv
    path "hcc1395_ml_summary.json",       emit: summary
    script:
    """
    python3 ${projectDir}/scripts/apply_model.py \\
        --vcf ${vcf} \\
        --modelos ${model_dir} \\
        --outdir .
    """
}
