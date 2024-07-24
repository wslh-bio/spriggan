process CALCULATE_ASSEMBLY {
    tag"$meta.id"
    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    label 'process_single'

    input:
    tuple val(meta), path(quast_report_tsv)
    path(kraken_results_tsv)
    path NCBI_assembly_stats_file

    output:
    path "*_Assembly_ratio_*"   , emit: assembly_ratio
    path "*_GC_content_*"   , emit: gc_content

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    calculate_assembly_ratio.py \\
    -d $NCBI_assembly_stats_file \\
    -q $quast_report_tsv \\
    -t $kraken_results_tsv
    """
}