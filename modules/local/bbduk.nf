process BBDUK {
    tag "$meta.id"
    label 'process_medium'

    container "staphb/bbtools:38.76"

    input:
    tuple val(meta), path(reads)
    path(contaminants)

    output:
    tuple val(meta), path('*trimmed*')      , emit: reads
    tuple val(meta), path('*.bbduk.log')    , emit: log
    path("*.trim.txt")                      , emit: bbduk_trim
    path("*.adapter.stats.txt")             , emit: bbduk_adapters
    path('*repaired*')                      , optional: true, emit: repaired_reads
    path('*_singletons.fastq.gz*')          , optional: true, emit: singletons
    path("*.repair.log")                    , optional: true, emit: repaired_log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def contaminants_fa = contaminants ? "ref=$contaminants" : "ref=/bbmap/resources/adapters.fa"
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')

    repair.sh \\
        -Xmx\$maxmem \\
        in1=${reads[0]} in2=${reads[1]} \\
        out1=${prefix}_repaired_1.fastq.gz out2=${prefix}_repaired_2.fastq.gz \\
        outs=${prefix}_singletons.fastq.gz \\
        repair \\
        &> ${prefix}.repair.log
        
    bbduk.sh \\
        -Xmx\$maxmem \\
        in1=${prefix}_repaired_1.fastq.gz in2=${prefix}_repaired_2.fastq.gz \\
        out1=${prefix}_trimmed_1.fastq.gz out2=${prefix}_trimmed_2.fastq.gz \\
        threads=$task.cpus \\
        $args \\
        $contaminants_fa \\
        stats=${prefix}.adapter.stats.txt \\
        &> ${prefix}.bbduk.log

    grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${prefix}.bbduk.log > ${prefix}.trim.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}