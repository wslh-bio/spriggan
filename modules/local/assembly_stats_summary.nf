process ASSEMBLY_STATS_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("data/*")

    output:
    path("assembly_stats_results.tsv"), emit: assembly_stats_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    assembly_stats_summary.py
    """
}
