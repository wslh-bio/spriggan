process COVERAGE_STATS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("data*/*")

    output:
    path('coverage_stats.tsv'), emit: coverage_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/python3.7

    import glob
    import os
    from numpy import median
    from numpy import average

    # function for summarizing samtools depth files
    def summarize_depth(file):
        # get sample id from file name and set up data list
        sid = os.path.basename(file).split('.')[0]
        data = []
        # open samtools depth file and get depth
        with open(file,'r') as inFile:
            for line in inFile:
                data.append(int(line.strip().split()[2]))
        # get median and average depth
        med = int(median(data))
        avg = int(average(data))
        # return sample id, median and average depth, and check for coverage fail
        if avg >= int(${params.mincoverage}):
            result = f"{sid}\\t{med}\\t{avg}\\tTRUE\\t\\n"
        if avg < int(${params.mincoverage}):
            result = f"{sid}\\t{med}\\t{avg}\\tFALSE\\tAverage coverage < ${params.mincoverage}X\\n"
        return result

    # get all samtools depth files
    files = glob.glob("data*/*.depth.tsv")

    # summarize samtools depth files
    results = map(summarize_depth,files)

    # write results to file
    with open('coverage_stats.tsv', 'w') as outFile:
        outFile.write("Sample\\tMedian Coverage\\tAverage Coverage\\tPass Coverage\\tComments\\n")
        for result in results:
            outFile.write(result)
    """
}
