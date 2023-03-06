process AMRFINDER_SETUP {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("kraken_results.tsv")
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.*.fa"), emit: amrfinder_input optional true

    script:
    """
    amrfinder_setup.py ${assembly}
    """
}
