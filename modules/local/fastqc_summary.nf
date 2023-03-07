process FASTQC_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/fastqc@sha256:6e4ace04f051e6138447ca64f5c6f44889157cc641fe8e0a40947a3949d47ab8"

    input:
    path('*')

    output:
    path("fastqc_summary.tsv"), emit: fastqc_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    zips=`ls *.zip`

    for i in \$zips; do
        unzip -o \$i &>/dev/null;
    done

    fq_folders=\${zips}

    for folder in \$fq_folders; do
        folder=\${folder%.*}
        cat \$folder/summary.txt >> fastqc_summary.tsv
        ls .
    done;

    sed -i 's/.fastq.gz//g' fastqc_summary.tsv
    """
}
