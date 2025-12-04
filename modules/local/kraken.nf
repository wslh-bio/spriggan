process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    container "staphb/kraken2:2.0.8-beta"

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("${meta.id}.kraken2.txt")   , emit: kraken_results
    path("kraken2.log")                               , optional: true, emit: log
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (params.kraken_db != null) {
        """
        dbname=${db}
        dbname=\${dbname%.*.*}

        mkdir custom-db
        tar -xvf ${db} --directory custom-db

        kraken2 --db ./custom-db --threads ${task.cpus} --report ${prefix}.kraken2.txt --paired ${reads[0]} ${reads[1]} | tee kraken2.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
            kraken DB: \$(echo \${db:-\$(echo \$dbname)})
        END_VERSIONS
        """
    } else {
        """
        kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${prefix}.kraken2.txt --paired ${reads[0]} ${reads[1]} | tee kraken2.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
            kraken DB: \$(echo \${db:-\$(ls /kraken2-db/)})
        END_VERSIONS
        """
    }
}
