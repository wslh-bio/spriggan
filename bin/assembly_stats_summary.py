#!/usr/bin/python3
import os
import glob
from pandas import DataFrame

# function for summarizing assembly output
def summarize_assembly_file(file):
    # get sample id from file name and set up data list
    sample_id = os.path.basename(file).split("_")[0]
    data = []
    data.append(sample_id)
    with open(file,"r") as inFile:
        for line in inFile:
            # Get Expected genome length
            if "Expected_length:" in line:
                expected_length = line.split(" ")[1].strip("\n")
                data.append(expected_length)
            # Get ratio
            if "Ratio Actual:Expected:" in line:
                ratio = line.split(" ")[2].strip("\n")
                data.append(ratio)
    return data

# get all calculate output files
assembly_files = glob.glob("data/*_Assembly_ratio_*")

# summarize output files
assembly_results = map(summarize_assembly_file, assembly_files)

# convert results to data frame and write to tsv
df = DataFrame(assembly_results,columns=['Sample','Expected Genome Length','Genome Length Ratio (Actual/Expected)'])

df.to_csv(f'assembly_stats_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')