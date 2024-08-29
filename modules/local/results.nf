process RESULTS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("bbduk_results.tsv")
    path("coverage_stats.tsv")
    path("quast_results.tsv")
    path("mlst_results.tsv")
    path("kraken_results.tsv")
    path("amrfinder_summary.tsv")
    path("selected_ar_genes.tsv")
    path(kraken_version, stageAs:"kraken_version.yml")
    path(amrfinder_version, stageAs:"amrfinder_version.yml")
    path("assembly_stats_summary.tsv")
    path("gc_stats_summary.tsv")

    output:
    path('spriggan_report.csv')

    script:
    """
    compile_results.py ${amrfinder_version} ${kraken_version} ${workflow.manifest.version}
    """
}
