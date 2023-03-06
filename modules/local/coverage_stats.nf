process COVERAGE_STATS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("data/*")

    output:
    path('coverage_stats.tsv'), emit: coverage_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    coverage_stats.py
    """
}
