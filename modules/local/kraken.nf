process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    container "staphb/kraken2:2.0.8-beta"

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.id}.kraken.txt")   , emit: kraken_results
    path("kraken.log")              , optional: true, emit: log
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kraken $args --db /kraken-database/minikraken_20171013_4GB --threads ${task.cpus} --paired ${reads[0]} ${reads[1]} > ${prefix}_raw.txt 2> kraken.log
    kraken-report --db /kraken-database/minikraken_20171013_4GB ${prefix}_raw.txt > ${prefix}.kraken.txt 2> kraken.log
    find -name kraken.log -size 0 -exec rm {} +

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken: \$(echo \$(kraken --version 2>&1) | sed 's/^.*Kraken //')
        kraken DB: \$(echo \$(ls /kraken-database/) )
    END_VERSIONS
    """
}
