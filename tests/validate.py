#!/usr/bin/env python3

import sys,os
import pandas as pd
import argparse

### Load validation data
ar_std = pd.read_csv("./validation/ar_summary_std.tsv",sep='\t',index_col="Sample")
cov_std = pd.read_csv("./validation/coverage_stats_std.tsv",sep='\t')
kraken_std = pd.read_csv("./validation/kraken_results_std.txt",sep='\t')
mlst_std = pd.read_csv("./validation/mlst_formatted_std.tsv",sep='\t')
quast_std = pd.read_csv("./validation/quast_results_std.tsv",sep='\t')

### Load in result data
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('ar_data',help='ar_summary.tsv')
parser.add_argument('cov_data',help='coverage_stats.tsv')
parser.add_argument('kraken_data',help='kraken_results.txt')
parser.add_argument('mlst_data',help='mlst_formatted.tsv')
parser.add_argument('quast_data',help='quast_results.tsv')
args = parser.parse_args()

ar_data = pd.read_csv(args.ar_data,sep='\t',index_col="Sample")

def check_coverage(y_std,x_data,range=2):
    x_data = round(float(x_data),2)
    y_std = round(float(y_std),2)
    if x_data >= y_std - 2 and x_data <= y_std + 2:
        return True
    else:
        return False

### Validate AR results
ar_hits = []
for sample in list(ar_std.index):

    if pd.isnull(ar_data.loc[sample].Gene):
        continue
    # get sample data and sort in list
    genes_data = ar_data.loc[sample].Gene.split(";")
    coverage_data = ar_data.loc[sample].Coverage.split(";")
    identity_data = ar_data.loc[sample].Identity.split(";")
    data = list(zip(genes_data,coverage_data,identity_data))
    data.sort()

    # get sample standard data and sort in list
    genes_std = ar_std.loc[sample].Gene.split(";")
    coverage_std = ar_std.loc[sample].Coverage.split(";")
    identity_std = ar_std.loc[sample].Identity.split(";")
    std = list(zip(genes_std,coverage_std,identity_std))
    std.sort()

    #compare data
    for a,b in zip(std,data):
        gene = a[0]
        coverage = check_coverage(a[1],b[1])
        identity = check_coverage(a[2],b[2])
        if not coverage:
            ar_hits.append({gene+"_coverage":" != ".join([a[1],b[1]])})
        if not identity:
            ar_hits.append({gene+"_identity":" != ".join([a[2],b[2]])})
