#!/usr/bin/python3
import os
import glob
import re

from pandas import DataFrame

# function for summarizing assembly output
def summarize_gc_file(file):
    # get sample id from file name and set up data list
    sample_id = os.path.basename(file).split("_")[0]
    data = []
    data.append(sample_id)
    with open(file,"r") as inFile:
        for line in inFile:
            # Get Species GC Content (mean)
            if "Species_GC_Mean:" in line:
                mean_GC_content = line.split(" ")[1].strip("\n")
                data.append(mean_GC_content)
            # Get GC percentage
            if "Sample_GC_Percent:" in line:
                sample_gc_content = line.split(" ")[1].strip("\n")
                data.append(sample_gc_content)
    return data

gc_files = glob.glob("data/*_GC_content_*")

gc_results = map(summarize_gc_file, gc_files)

df = DataFrame(gc_results,columns=['Sample', 'Species GC Content (Mean)', 'Sample GC Content (%)'])

df.to_csv(f'gc_stats_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')