process MUTECT2 {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    label 'process_gatk'
    conda 'bioconda::gatk4=4.5.0.0'
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(tumor_meta),  path(tumor_bam),  path(tumor_bai)
    tuple val(normal_meta), path(normal_bam), path(normal_bai)
    path genome
    path genome_index

    output:
    path "somatic.vcf.gz",       emit: vcf
    path "somatic.vcf.gz.tbi",   emit: tbi
    path "somatic.vcf.gz.stats", emit: stats

    script:
    """
    gatk Mutect2 \\
        -R ${genome} \\
        -I ${tumor_bam}  --tumor-sample  ${tumor_meta.id} \\
        -I ${normal_bam} --normal-sample ${normal_meta.id} \\
        -O somatic.vcf.gz

    gatk FilterMutectCalls \\
        -R ${genome} \\
        -V somatic.vcf.gz \\
        --stats somatic.vcf.gz.stats \\
        -O somatic.filtered.vcf.gz
    """
}
