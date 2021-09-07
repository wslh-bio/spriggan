#!/usr/bin/env python3

import sys,os
import pandas as pd
import argparse

base_path =  os.path.dirname(os.path.realpath(__file__))

### Load validation data
ar_std = pd.read_csv(os.path.join(base_path,"validation/ar_summary_std.tsv"),sep='\t',index_col="Sample")
cov_std = pd.read_csv(os.path.join(base_path,"validation/coverage_stats_std.tsv"),sep='\t',index_col="Sample")
kraken_std = pd.read_csv(os.path.join(base_path,"validation/kraken_results_std.tsv"),sep='\t',index_col="Sample")
mlst_std = pd.read_csv(os.path.join(base_path,"validation/mlst_results_std.tsv"),sep='\t',index_col="Sample")
quast_std = pd.read_csv(os.path.join(base_path,"validation/quast_results_std.tsv"),sep='\t',index_col="Sample")

### Load in result data
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('ar_data',help='ar_summary.tsv')
parser.add_argument('cov_data',help='coverage_stats.tsv')
parser.add_argument('kraken_data',help='kraken_results.txt')
parser.add_argument('mlst_data',help='mlst_formatted.tsv')
parser.add_argument('quast_data',help='quast_results.tsv')
args = parser.parse_args()

ar_data = pd.read_csv(args.ar_data,sep='\t',index_col="Sample")
cov_data = pd.read_csv(args.cov_data,sep='\t',index_col="Sample")
kraken_data = pd.read_csv(args.kraken_data,sep='\t',index_col="Sample")
mlst_data = pd.read_csv(args.mlst_data,sep='\t',index_col="Sample")
quast_data = pd.read_csv(args.quast_data,sep='\t',index_col="Sample")

def check_compare(y_std,x_data,range=2):
    x_data = round(float(x_data),2)
    y_std = round(float(y_std),2)
    if x_data >= y_std - range and x_data <= y_std + range:
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
        coverage = check_compare(a[1],b[1])
        identity = check_compare(a[2],b[2])
        if not coverage:
            ar_hits.append({sample+"_"+gene+"_coverage":" != ".join([a[1],b[1]])})
        if not identity:
            ar_hits.append({sample+"_"+gene+"_identity":" != ".join([a[2],b[2]])})

### Validate Coverage results
cov_hits = []
for sample in list(cov_std.index):
    med_cov = check_compare(cov_std.loc[sample]["Median Coverage"],cov_data.loc[sample]["Median Coverage"])
    avg_cov = check_compare(cov_std.loc[sample]["Average Coverage"],cov_data.loc[sample]["Average Coverage"])
    if not med_cov:
        cov_hits.append({sample+"_"+"Median Coverage":" != ".join([str(cov_std.loc[sample]["Median Coverage"]),str(cov_data.loc[sample]["Median Coverage"])])})
    if not avg_cov:
        cov_hits.append({sample+"_"+"Average Coverage":" != ".join([str(cov_std.loc[sample]["Average Coverage"]),str(cov_data.loc[sample]["Average Coverage"])])})

### Validate Kraken results
kraken_hits = []
for sample in list(kraken_std.index):
    unclass_pct_std = kraken_std.loc[sample]["Unclassified Reads (%)"].strip("%")
    primary_name_std = kraken_std.loc[sample]["Primary Species (%)"].split('(')[0]
    primary_pct_std = kraken_std.loc[sample]["Primary Species (%)"].split('(')[1].strip("()%")
    secondary_name_std = kraken_std.loc[sample]["Secondary Species (%)"].split('(')[0]
    secondary_pct_std = kraken_std.loc[sample]["Secondary Species (%)"].split('(')[1].strip("()%")

    unclass_pct = kraken_data.loc[sample]["Unclassified Reads (%)"].strip("%")
    primary_name = kraken_data.loc[sample]["Primary Species (%)"].split('(')[0]
    primary_pct = kraken_data.loc[sample]["Primary Species (%)"].split('(')[1].strip("()%")
    secondary_name = kraken_data.loc[sample]["Secondary Species (%)"].split('(')[0]
    secondary_pct = kraken_data.loc[sample]["Secondary Species (%)"].split('(')[1].strip("()%")

    primary_pct_result = check_compare(primary_pct_std,primary_pct)
    secondary_pct_result = check_compare(secondary_pct_std,secondary_pct)
    unclass_pct_result = check_compare(unclass_pct_std,unclass_pct)

    if not primary_pct_result:
        kraken_hits.append({sample+"-Primary":" != ".join([primary_name_std+" @ "+primary_pct_std,primary_name+" @ "+primary_pct])})
    if not secondary_pct_result:
        kraken_hits.append({sample+"-Secondary":" != ".join([secondary_name_std+" @ "+secondary_pct_std,secondary_name+" @ "+secondary_pct])})
    if not unclass_pct_result:
        kraken_hits.append({sample+"-Unclassified":" != ".join([unclass_pct_std,unclass_pct])})

### Validate MLST results
mlst_hits = []
for sample in list(mlst_std.index):
    if ";" in mlst_data.loc[sample]["MLST Scheme"]:
        std = mlst_std.loc[sample]["MLST Scheme"].split(";")
        test = mlst_data.loc[sample]["MLST Scheme"].split(";")
        std.sort()
        test.sort()
        std = ";".join(std)
        test = ";".join(test)
    else:
        std = mlst_std.loc[sample]["MLST Scheme"]
        test = mlst_data.loc[sample]["MLST Scheme"]

    if std != test:
        mlst_hits.append({sample:" != ".join([std,test])})

### Validate QUAST results
quast_hits = []
for sample in list(quast_std.index):

    contigs_std = str(quast_std.loc[sample]["Contigs"])
    length_std = str(quast_std.loc[sample]["Assembly Length (bp)"])
    n50_std = str(quast_std.loc[sample]["N50"])

    contigs = str(quast_data.loc[sample]["Contigs"])
    length = str(quast_data.loc[sample]["Assembly Length (bp)"])
    n50 = str(quast_data.loc[sample]["N50"])

    contig_result = check_compare(contigs_std,contigs,10)
    length_result = check_compare(length_std,length,50000)
    n50_result = check_compare(n50_std,n50,500)

    if not contig_result:
        quast_hits.append({sample+"_contigs":contigs_std+" != "+contigs})
    if not length_result:
        quast_hits.append({sample+"_length":length_std+" != "+length})
    if not n50_result:
        quast_hits.append({sample+"_N50":n50_std+" != "+n50})

### Report Results
validation_pass = True

# AR
print("AR Validation")
if not ar_hits:
    print("--Pass--")
else:
    for hit in ar_hits:
        print(hit)
    validation_pass = False

# Coverage
print("Coverage Validation")
if not cov_hits:
    print("--Pass--")
else:
    for hit in cov_hits:
        print(hit)
    validation_pass = False

# Kraken
print("Kraken Validation")
if not kraken_hits:
    print("--Pass--")
else:
    for hit in kraken_hits:
        print(hit)
    validation_pass = False

# MLST
print("MLST Validation")
if not mlst_hits:
    print("--Pass--")
else:
    for hit in mlst_hits:
        print(hit)
    validation_pass = False

# QUAST
print("QUAST Validation")
if not quast_hits:
    print("--Pass--")
else:
    for hit in quast_hits:
        print(hit)
    validation_pass = False

if not validation_pass:
    sys.exit(1)
