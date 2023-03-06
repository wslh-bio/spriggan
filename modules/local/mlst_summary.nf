process MLST_SUMMARY {
    label 'process_single'


    container 'quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2'

    input:
    path("data/*")

    output:
    path("mlst_results.tsv"), emit: mlst_tsv

    script:
    """
    mlst_summary.py
    """
}
