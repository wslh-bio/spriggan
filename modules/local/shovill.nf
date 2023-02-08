process SHOVILL {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    container "staphb/shovill:1.1.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa")                              , emit: contigs
    tuple val(meta), path("${meta.id}.sam")                                     , emit: sam_files
    tuple val(meta), path("shovill_output/shovill.corrections")                 , emit: corrections
    tuple val(meta), path("shovill_output/shovill.log")                         , emit: log
    tuple val(meta), path("shovill_output/{skesa,spades,megahit,velvet}.fasta") , emit: raw_contigs
    tuple val(meta), path("shovill_output/contigs.{fastg,gfa,LastGraph}")       , optional:true, emit: gfa
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory.toGiga()
    """
    shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        $args \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./shovill_output \\
        --force

    mv shovill_output/contigs.fa ${prefix}.contigs.fa
    bwa index ${prefix}.contigs.fa
    bwa mem ${prefix}.contigs.fa ${reads[0]} ${reads[1]} > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(echo \$(shovill --version 2>&1) | sed 's/^.*shovill //')
    END_VERSIONS
    """
}
