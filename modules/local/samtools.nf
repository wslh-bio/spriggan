process SAMTOOLS {
    tag "$meta.id"
    label 'process_medium'

    container "staphb/samtools:1.10"

    input:
    tuple val(meta), path(sam_files)

    output:
    path("*.depth.tsv")     , emit: cov_files
    path("*.stats.txt")     , emit: stats_multiqc
    path("*.bam")           , emit: sorted_bam
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory.toGiga()
    """
    samtools view -S -b ${prefix}.sam | samtools sort > ${prefix}.bam
    samtools index ${prefix}.bam
    samtools depth -a ${prefix}.bam > ${prefix}.depth.tsv
    samtools stats ${prefix}.bam > ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //')
    END_VERSIONS
    """
}
