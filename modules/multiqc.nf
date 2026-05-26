process MULTIQC {
    tag "multiqc"
    label 'low_cpu'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container 'multiqc/multiqc:v1.21'
    input:
    path zips
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
    script:
    """
    multiqc . --filename multiqc_report.html
    """
}
