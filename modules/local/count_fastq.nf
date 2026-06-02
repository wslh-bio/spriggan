process COUNT_FASTQ {
    tag "$meta.id"
    label 'process_single'
    container "quay.io/wslh-bioinformatics/pandas:1.5.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    count_fastq.py \\
        $args \\
        --sample_id ${meta.id} \\
        --raw_reads ${reads} \\
        --output ${meta.id}_output.csv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        countfastq: \$(python --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}_output.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        countfastq: \$(python --version)
    END_VERSIONS
    """
}
