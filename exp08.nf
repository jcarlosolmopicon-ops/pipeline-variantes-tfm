#!/usr/bin/env nextflow
/*
 * exp08.nf — Inferencia del modulo de clasificacion como proceso de Nextflow.
 *
 * Workflow independiente de main.nf: aplica los modelos serializados de
 * exp04 al VCF anotado por VEP mediante el proceso CLASSIFY, sin repetir
 * ninguna etapa del pipeline de deteccion.
 *
 * Uso:
 *   nextflow run exp08.nf -profile conda_local -c config/exp08.config \
 *     --vcf resultados/exp03/vep/annotated.vcf.gz \
 *     --modelos resultados/exp04 \
 *     --outdir resultados/exp08
 */

nextflow.enable.dsl = 2

include { CLASSIFY } from './modules/classify'

params.vcf     = null
params.modelos = null

workflow {
    if (!params.vcf || !params.modelos) {
        error "Faltan parametros: --vcf <annotated.vcf.gz> --modelos <dir con model_*.joblib>"
    }

    log.info """
╔══════════════════════════════════════════════════════════════╗
║   exp08 — CLASIFICACION ML COMO PROCESO DE NEXTFLOW          ║
╚══════════════════════════════════════════════════════════════╝
 vcf      : ${params.vcf}
 modelos  : ${params.modelos}
 outdir   : ${params.outdir}
 profile  : ${workflow.profile}
──────────────────────────────────────────────────────────────
""".stripIndent()

    vcf_ch    = Channel.value(file(params.vcf, checkIfExists: true))
    models_ch = Channel.value(file(params.modelos, checkIfExists: true))

    CLASSIFY(vcf_ch, models_ch)
}
