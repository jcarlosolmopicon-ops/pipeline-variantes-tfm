#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC      } from './modules/fastqc'
include { MULTIQC     } from './modules/multiqc'
include { TRIM_GALORE } from './modules/trim_galore'
include { BWA_MEM2    } from './modules/bwa_mem2'
include { MARKDUPLICATES } from './modules/markduplicates'
include { BQSR } from './modules/bqsr'
include { MUTECT2 } from './modules/mutect2'
include { VEP     } from './modules/vep'

workflow {
    log.info """
╔══════════════════════════════════════════════════════════════╗
║       PIPELINE VARIANTES SOMÁTICAS — TFM UAX 2026           ║
╚══════════════════════════════════════════════════════════════╝
 input    : ${params.input}
 outdir   : ${params.outdir}
 genome   : ${params.genome}
 profile  : ${workflow.profile}
 step     : ${params.step}
──────────────────────────────────────────────────────────────
""".stripIndent()

    Channel
        .fromFilePairs(params.input, checkIfExists: true)
        .map { id, reads -> [ [id: id], reads ] }
        .set { reads_ch }

    if (params.step == 'qc' || params.step == 'all') {
        FASTQC(reads_ch)
        MULTIQC(FASTQC.out.zip.map { it[1] }.collect())
    }

    if (params.step == 'trim' || params.step == 'all') {
        TRIM_GALORE(reads_ch)
    }

    if (params.step == 'align' || params.step == 'all') {
        genome_ch  = Channel.value(file(params.genome))
        index_ch   = Channel.fromPath("${params.genome}.*").collect()
        dict_ch    = Channel.value(file(params.genome.replace('.fa', '.dict')))
        dbsnp_ch   = Channel.value(file(params.dbsnp))
        dbsnp_idx  = Channel.value(file(params.dbsnp_index))
        TRIM_GALORE(reads_ch)
        BWA_MEM2(TRIM_GALORE.out.trimmed_reads, genome_ch, index_ch)
        MARKDUPLICATES(BWA_MEM2.out.bam.join(BWA_MEM2.out.bai))
        BQSR(
            MARKDUPLICATES.out.bam.join(MARKDUPLICATES.out.bai),
            genome_ch,
            index_ch,
            dict_ch,
            dbsnp_ch,
            dbsnp_idx
        )
    }

    if (params.step == 'call' || params.step == 'all') {
        genome_ch  = Channel.value(file(params.genome))
        index_ch   = Channel.fromPath("${params.genome}.*").collect()
        dict_ch    = Channel.value(file(params.genome.replace('.fa', '.dict')))

        tumor_ch  = Channel.fromPath("${params.outdir}/bqsr/tumor.bqsr.bam")
            .map { bam -> [ [id: 'tumor'],  bam, file("${bam}.bai") ] }
        normal_ch = Channel.fromPath("${params.outdir}/bqsr/normal.bqsr.bam")
            .map { bam -> [ [id: 'normal'], bam, file("${bam}.bai") ] }

        MUTECT2(tumor_ch, normal_ch, genome_ch, index_ch, dict_ch)
    }

    if (params.step == 'annotate' || params.step == 'all') {
        vcf_ch = Channel.fromPath("${params.outdir}/mutect2/somatic.filtered.vcf.gz")
        VEP(vcf_ch)
    }
}
