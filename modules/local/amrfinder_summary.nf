process AMRFINDER_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("data/*")

    output:
    path("amrfinder_predictions.tsv")
    path("amrfinder_summary.tsv"), emit: amrfinder_tsv
    path("selected_ar_genes.tsv"), emit: selected_ar_tsv

    script:
    """
    amrfinder_summary.py
    """
}
