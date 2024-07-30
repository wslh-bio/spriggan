process GC_STATS_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("data/*")

    output:
    path("gc_stats_results.tsv"), emit: gc_stats_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gc_stats_summary.py
    """
}
