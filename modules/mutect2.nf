process MUTECT2 {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    label 'gatk'
    publishDir "${params.outdir}/mutect2", mode: 'copy'
    conda 'bioconda::gatk4=4.5.0.0'
    errorStrategy 'ignore'
    input:
    tuple val(tumor_meta),  path(tumor_bam),  path(tumor_bai)
    tuple val(normal_meta), path(normal_bam), path(normal_bai)
    path genome
    path genome_index
    path genome_dict
    output:
    path "somatic.vcf.gz",          emit: vcf
    path "somatic.vcf.gz.tbi",      emit: tbi
    path "somatic.vcf.gz.stats",    emit: stats
    path "somatic.filtered.vcf.gz", emit: vcf_filtered
    script:
    """
    gatk --java-options "-Xmx24g" Mutect2 \\
        -R ${genome} \\
        -I ${tumor_bam}  --tumor-sample  ${tumor_meta.id} \\
        -I ${normal_bam} --normal-sample ${normal_meta.id} \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O somatic.vcf.gz
    gatk --java-options "-Xmx24g" FilterMutectCalls \\
        -R ${genome} \\
        -V somatic.vcf.gz \\
        --stats somatic.vcf.gz.stats \\
        -O somatic.filtered.vcf.gz
    """
}
