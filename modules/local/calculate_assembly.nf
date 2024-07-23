process CALCULATE_ASSEMBLY {
    tag "$each.id"
    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    tuple val(each), path(quast_report_tsv)
    path NCBI_assembly_stats_file
    path kraken_results_tsv

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
    -t $kraken_results_tsv \\
    -s $each.id
    """
}