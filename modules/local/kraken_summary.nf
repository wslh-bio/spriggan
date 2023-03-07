process KRAKEN_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("data/*")

    output:
    path("kraken_results.tsv"), emit: kraken_tsv

    script:
    """
    kraken_summary.py
    """
}
