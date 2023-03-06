#!/usr/bin/python3
import os
import glob
from pandas import DataFrame

# function for summarizing bbduk output
def summarize_bbduk(file):
    # get sample id from file name and set up data list
    sample_id = os.path.basename(file).split(".")[0]
    data = []
    data.append(sample_id)
    with open(file,"r") as inFile:
        for line in inFile:
            # get total number of reads
            if "Result:" in line:
                num_reads = line.split("reads")[0].strip()
                data.append(num_reads)
            # get total number of reads removed
            if "Total Removed:" in line:
                rm_reads = line.split("reads")[0].strip()
                rm_reads = rm_reads.rstrip()
                data.append(rm_reads)
    return data

# get all bbduk output files
files = glob.glob("data/*.trim.txt")

# summarize bbduk output files
results = map(summarize_bbduk,files)

# convert results to data frame and write to tsv
df = DataFrame(results,columns=['Sample','Reads Removed','Total Reads'])
print(df)
df.to_csv(f'bbduk_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')