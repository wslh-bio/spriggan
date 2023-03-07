process BBDUK_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("data/*")

    output:
    path("bbduk_results.tsv"), emit: bbduk_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bbduk_summary.py
    """
}
