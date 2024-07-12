process CALCULATE_RATIO {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    path fastas
    path reference_fasta


    output:
    path ""   , emit: 

    script: // This script is bundled with the pipeline, in wslh-bio/dryad/bin/
    """
    calculate_ratio.py \\
    $fastas 
    """
}