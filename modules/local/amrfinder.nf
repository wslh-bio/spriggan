process AMRFINDER {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/wslh-bioinformatics/spriggan-amrfinder:3.10.30"

    input:
    tuple val(meta), path(amrfinder_input)

    output:
    path("${meta.id}.amr.tsv"), emit: amrfinder_predictions
    path "versions.yml"    , emit: versions

    script:
    """
    run_amrfinder.py ${amrfinder_input} ${params.plus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinder: \$(echo \$(amrfinder --version) )
        amrfinder DB: \$(echo \$(cat /amrfinder/data/latest/version.txt) )
    END_VERSIONS
    """
}
