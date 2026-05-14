#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC      } from './modules/fastqc'
include { MULTIQC     } from './modules/multiqc'
include { TRIM_GALORE } from './modules/trim_galore'
include { BWA_MEM2    } from './modules/bwa_mem2'

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
        genome_ch = Channel.value(file(params.genome))
        index_ch  = Channel.fromPath("${params.genome}.*").collect()
        TRIM_GALORE(reads_ch)
        BWA_MEM2(TRIM_GALORE.out.trimmed_reads, genome_ch, index_ch)
    }
}
