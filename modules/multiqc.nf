process MULTIQC {
    label 'process_low'

    conda 'bioconda::multiqc=1.21'
    container 'multiqc/multiqc:v1.21'

    input:
    path(reports)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data/"),       emit: data

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}
