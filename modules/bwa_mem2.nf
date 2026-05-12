process BWA_MEM2 {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19'
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    tuple val(meta), path(reads)
    path genome
    path genome_index

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.sorted.bam.bai"), emit: bai

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}_lib1"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${genome} \\
        ${reads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -o ${meta.id}.sorted.bam

    samtools index ${meta.id}.sorted.bam
    """
}
