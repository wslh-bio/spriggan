process QUAST {
    tag "$meta.id"
    label 'process_medium'

    container "staphb/quast:5.0.2"

    input:
    tuple val(meta), path(contigs)
    val(min_quast_contig)

    output:
    tuple val(meta), path("${meta.id}.transposed.quast.report.tsv") , emit: transposed_report
    tuple val(meta), path("${meta.id}.quast.report.tsv")            , emit: result
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    quast.py ${contigs} -m ${min_quast_contig} -o .
    mv report.tsv ${prefix}.quast.report.tsv
    mv transposed_report.tsv ${prefix}.transposed.quast.report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
