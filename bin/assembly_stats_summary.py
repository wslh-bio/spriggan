#!/usr/bin/python3
import os
import glob
from pandas import DataFrame
import sys

# function for summarizing bbduk output
def summarize_assembly_file(file):
    # get sample id from file name and set up data list
    sample_id = os.path.basename(file).split("_")[0]
    data = []
    data.append(sample_id)
    with open(file,"r") as inFile:
        for line in inFile:
            # get total number of reads
            if "Ratio Actual:Expected:" in line:
                ratio = line.split(" ")[2].strip("\n")
                data.append(ratio)
    return data

# get all bbduk output files
files = glob.glob("data/*_Assembly_ratio_*")

# summarize bbduk output files
results = map(summarize_assembly_file,files)

# convert results to data frame and write to tsv
df = DataFrame(results,columns=['Sample','Ratio Actual:Expected'])
print(df)
df.to_csv(f'assembly_stats_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')