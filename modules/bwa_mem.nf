process BWA_MEM {
    tag "${meta.id}"
    label 'high_cpu'
    publishDir "${params.outdir}/bam", mode: 'copy'
    container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'
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
    bwa mem -t ${task.cpus} -R "${rg}" ${genome} ${reads} | samtools sort -@ ${task.cpus} -m 1500M -o ${meta.id}.sorted.bam
    samtools index ${meta.id}.sorted.bam
    """
}
