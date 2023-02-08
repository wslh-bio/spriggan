process BBDUK_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("data*/*")

    output:
    path("bbduk_results.tsv"), emit: bbduk_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/python3
    import os
    import glob
    import numpy
    import pandas as pd
    from pandas import DataFrame

    # function for summarizing bbduk output
    def summarize_bbduk(file):
        # get sample id from file name and set up data list
        sample_id = os.path.basename(file).split(".")[0]
        data = []
        data.append(sample_id)
        with open(file,"r") as inFile:
            for i, line in enumerate(inFile):
                # get total number of reads
                if i == 0:
                    num_reads = line.strip().split("\\t")[1].replace(" reads ","")
                    data.append(num_reads)
                # get total number of reads removed
                if i == 3:
                    rm_reads = line.strip().split("\\t")[1].replace("reads ","")
                    rm_reads = rm_reads.rstrip()
                    data.append(rm_reads)
        return data

    # get all bbduk output files
    files = glob.glob("data*/*.trim.txt")

    # summarize bbduk output files
    results = map(summarize_bbduk,files)

    # convert results to data frame and write to tsv
    df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])
    df.to_csv(f'bbduk_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
    """
}
