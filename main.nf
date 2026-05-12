#!/usr/bin/env nextflow

/*
 * main.nf — Pipeline de variantes somáticas (tumor vs normal)
 * TFM: Juan Carlos Olmo Picón — UAX Bioinformática 2025
 *
 * Dataset:    SEQC2/HCC1395  (exp01 — validación)
 *             TCGA-LUAD MAF  (exp02 — módulo ML)
 * Referencia: GRCh38/hg38
 *
 * Uso:
 *   nextflow run main.nf -profile conda_local
 *   nextflow run main.nf -profile test
 *   nextflow run main.nf -profile conda_local -entry QC
 *
 * Refs: NB-005, NB-006, NXF-001, ADR-01, ADR-03, ADR-04, ADR-05
 */

nextflow.enable.dsl = 2

log.info """
╔══════════════════════════════════════════════════════════════╗
║       PIPELINE VARIANTES SOMÁTICAS — TFM UAX 2026           ║
╚══════════════════════════════════════════════════════════════╝
 input    : ${params.input}
 outdir   : ${params.outdir}
 genome   : ${params.genome}
 profile  : ${workflow.profile}
 version  : ${manifest.version}
──────────────────────────────────────────────────────────────
""".stripIndent()


// ─── MÓDULO: FastQC ───────────────────────────────────────────────────────────

process FASTQC {
    tag         "$sample_id"
    label       'low_cpu'
    publishDir  "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads}
    """

    stub:
    """
    touch ${sample_id}_R1_fastqc.html ${sample_id}_R1_fastqc.zip
    touch ${sample_id}_R2_fastqc.html ${sample_id}_R2_fastqc.zip
    """
}


// ─── MÓDULO: MultiQC ──────────────────────────────────────────────────────────

process MULTIQC {
    label       'low_cpu'
    publishDir  "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(fastqc_zips)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc \\
        --title "TFM Variantes Somáticas — Control de Calidad" \\
        --outdir . \\
        .
    """

    stub:
    """
    mkdir -p multiqc_data
    touch multiqc_report.html multiqc_data/multiqc_general_stats.txt
    """
}


// ─── MÓDULO: Trim Galore ──────────────────────────────────────────────────────

process TRIM_GALORE {
    tag         "$sample_id"
    label       'low_cpu'
    publishDir  "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_{1,2}.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt",                      emit: trim_log

    script:
    """
    trim_galore \\
        --paired \\
        --cores ${task.cpus} \\
        --gzip \\
        ${reads}
    """

    stub:
    """
    touch ${sample_id}_val_1.fq.gz ${sample_id}_val_2.fq.gz
    touch ${sample_id}_trimming_report.txt
    """
}


// ─── MÓDULO: BWA-MEM2 + SAMtools sort/index ───────────────────────────────────

process BWA_MEM2 {
    tag         "$sample_id"
    label       'high_cpu'
    publishDir  "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path genome
    path genome_index

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"),     emit: bam
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai

    script:
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}_lib1"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${genome} \\
        ${reads} \\
    | samtools sort \\
        -@ ${task.cpus} \\
        -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
    """
}


// ─── MÓDULO: GATK MarkDuplicates ──────────────────────────────────────────────

process MARKDUPLICATES {
    tag         "$sample_id"
    label       'gatk'
    publishDir  "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"),     emit: bam
    tuple val(sample_id), path("${sample_id}.markdup.bam.bai"), emit: bai
    path "${sample_id}.markdup.metrics",                         emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.markdup.bam \\
        -M ${sample_id}.markdup.metrics \\
        --TMP_DIR /tmp

    samtools index ${sample_id}.markdup.bam
    """

    stub:
    """
    touch ${sample_id}.markdup.bam ${sample_id}.markdup.bam.bai ${sample_id}.markdup.metrics
    """
}


// ─── MÓDULO: GATK BQSR ────────────────────────────────────────────────────────

process BQSR {
    tag         "$sample_id"
    label       'gatk'
    publishDir  "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome
    path genome_index

    output:
    tuple val(sample_id), path("${sample_id}.bqsr.bam"),     emit: bam
    tuple val(sample_id), path("${sample_id}.bqsr.bam.bai"), emit: bai

    script:
    """
    gatk BaseRecalibrator \\
        -I ${bam} \\
        -R ${genome} \\
        --known-sites /mnt/f/TFM-bioinformatica/datos-raw/referencia/dbsnp_156.vcf.gz \\
        -O ${sample_id}.recal.table

    gatk ApplyBQSR \\
        -I ${bam} \\
        -R ${genome} \\
        --bqsr-recal-file ${sample_id}.recal.table \\
        -O ${sample_id}.bqsr.bam

    samtools index ${sample_id}.bqsr.bam
    """

    stub:
    """
    touch ${sample_id}.bqsr.bam ${sample_id}.bqsr.bam.bai
    """
}


// ─── MÓDULO: GATK Mutect2 (tumor vs normal) ───────────────────────────────────

process MUTECT2 {
    tag         "tumor_vs_normal"
    label       'gatk'
    publishDir  "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(tumor_id),  path(tumor_bam),  path(tumor_bai)
    tuple val(normal_id), path(normal_bam), path(normal_bai)
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
        -I ${tumor_bam}  --tumor-sample  ${tumor_id} \\
        -I ${normal_bam} --normal-sample ${normal_id} \\
        -O somatic.vcf.gz

    gatk FilterMutectCalls \\
        -R ${genome} \\
        -V somatic.vcf.gz \\
        --stats somatic.vcf.gz.stats \\
        -O somatic.filtered.vcf.gz
    """

    stub:
    """
    touch somatic.vcf.gz somatic.vcf.gz.tbi somatic.vcf.gz.stats
    """
}


// ─── MÓDULO: VEP ──────────────────────────────────────────────────────────────

process VEP {
    tag         "vep_annotation"
    label       'vep'
    publishDir  "${params.outdir}/vep", mode: 'copy'

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

    stub:
    """
    touch annotated.vcf.gz annotated.vcf.gz_summary.html
    """
}


// ─── ENTRADA: Solo QC (NXF-001) ───────────────────────────────────────────────

workflow QC {
    Channel
        .fromFilePairs(params.input, checkIfExists: true)
        .set { reads_ch }

    FASTQC(reads_ch)
    MULTIQC(FASTQC.out.zip.map { it[1] }.collect())

    emit:
    fastqc_html = FASTQC.out.html
    multiqc     = MULTIQC.out.report
}


// ─── WORKFLOW PRINCIPAL ───────────────────────────────────────────────────────

workflow {

    Channel
        .fromFilePairs(params.input, checkIfExists: true)
        .set { raw_reads_ch }

    FASTQC(raw_reads_ch)
    MULTIQC(FASTQC.out.zip.map { it[1] }.collect())

    TRIM_GALORE(raw_reads_ch)

    genome_ch = Channel.value(file(params.genome))
    index_ch  = Channel.value(file(params.genome_index))

    BWA_MEM2(TRIM_GALORE.out.trimmed_reads, genome_ch, index_ch)

    MARKDUPLICATES(BWA_MEM2.out.bam.join(BWA_MEM2.out.bai))
    BQSR(
        MARKDUPLICATES.out.bam.join(MARKDUPLICATES.out.bai),
        genome_ch,
        index_ch
    )

    // Mutect2 y VEP se activan en NXF-002/003
    // una vez que los FASTQs de HCC1395 estén convertidos
    // MUTECT2(tumor_ch, normal_ch, genome_ch, index_ch)
    // VEP(MUTECT2.out.vcf)
}
