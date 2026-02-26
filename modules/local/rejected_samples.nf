process REJECTED_SAMPLES {

    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    path csv_file
    val workflow_name

    output:
    path "*_empty_samples.csv", emit: rejected_samples_report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash
    cp ${csv_file} ${workflow_name}_empty_samples.csv
    """
}