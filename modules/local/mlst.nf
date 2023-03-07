process MLST {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/wslh-bioinformatics/spriggan-mlst:2.19.0"

    input:
    tuple val(meta), path(input)

    output:
    path("*.mlst.tsv")     , emit: mlst_files
    path "versions.yml"    , emit: versions
    path("*.alleles.tsv")


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_mlst.py ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
    END_VERSIONS
    """
}
