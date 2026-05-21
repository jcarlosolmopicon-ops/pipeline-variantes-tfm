process BWA_MEM2 {
    tag "${meta.id}"
    label 'high_cpu'
    publishDir "${params.outdir}/bam", mode: 'copy'
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.19'
    errorStrategy 'ignore'
    input:
    tuple val(meta), path(reads)
    path genome
    path index
    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai
    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}_lib1"
    """
    bwa mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${genome} \\
        ${reads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -m 4G \\
        -o ${meta.id}.sorted.bam
    samtools index ${meta.id}.sorted.bam
    """
}
