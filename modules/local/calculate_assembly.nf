process CALCULATE_ASSEMBLY {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    path NCBI_assembly_stats_file
    path quast_report_tsv
    path kraken_results_tsv
    val sample



    output:
    path "*_Assembly_ratio_*"   , emit: assembly_ratio
    path "*GC_content_*"   , emit: gc_content

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    calculate_assembly_ratio.py \\
    -d $NCBI_assembly_stats_file \\
    -q $quast_report_tsv \\
    -t $kraken_results_tsv \\
    -s $sample
    """
}