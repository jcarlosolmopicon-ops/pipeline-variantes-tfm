#!/usr/bin/env nextflow

/*
 * main.nf — Pipeline de variantes somáticas (tumor vs normal)
 * TFM: Juan Carlos Olmo Picón — UAX Bioinformática 2026
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

include { FASTQC        } from './modules/fastqc'
include { MULTIQC       } from './modules/multiqc'
include { TRIM_GALORE   } from './modules/trim_galore'
include { BWA_MEM2      } from './modules/bwa_mem2'
include { MARKDUPLICATES} from './modules/markduplicates'
include { BQSR          } from './modules/bqsr'
include { MUTECT2       } from './modules/mutect2'
include { VEP           } from './modules/vep'

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


// ─── ENTRADA: Solo QC (NXF-001) ───────────────────────────────────────────────

workflow QC {
    Channel
        .fromFilePairs(params.input, checkIfExists: true)
        .map { id, reads -> [ [id: id], reads ] }
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
        .map { id, reads -> [ [id: id], reads ] }
        .set { raw_reads_ch }

    // 1. QC inicial
    FASTQC(raw_reads_ch)
    MULTIQC(FASTQC.out.zip.map { it[1] }.collect())

    // 2. Trimado
    TRIM_GALORE(raw_reads_ch)

    // 3. Alineamiento
    genome_ch = Channel.value(file(params.genome))
    index_ch  = Channel.value(file(params.genome_index))

    BWA_MEM2(TRIM_GALORE.out.trimmed_reads, genome_ch, index_ch)

    // 4. MarkDuplicates + BQSR
    MARKDUPLICATES(BWA_MEM2.out.bam.join(BWA_MEM2.out.bai))

    dbsnp_ch = Channel.value(file("/mnt/f/TFM-bioinformatica/datos-raw/referencia/dbsnp_156.vcf.gz"))

    BQSR(
        MARKDUPLICATES.out.bam.join(MARKDUPLICATES.out.bai),
        genome_ch,
        index_ch,
        dbsnp_ch
    )

    // 5. Mutect2 y VEP — activar en NXF-002/003
    //    una vez que los FASTQs de HCC1395 estén convertidos
    // MUTECT2(tumor_ch, normal_ch, genome_ch, index_ch)
    // VEP(MUTECT2.out.vcf)
}
