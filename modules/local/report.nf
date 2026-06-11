process REPORT {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path results_compiled
    val empty_ntc

    output:
    path("${params.run_name}_spriggan_report.csv"), emit: result_csv

    script:
    """
    create_report.py \
        --result_files ${results_compiled} \
        --workflowVersion ${workflow.manifest.version} \
        --workflowRunName ${params.run_name} \
        --empty_ntc_list ${empty_ntc}
    """
}
