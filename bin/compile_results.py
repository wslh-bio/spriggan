#!/usr/bin/env python3
import sys
import glob
import pandas as pd
from functools import reduce

with open(sys.argv[1], 'r') as amrFile:
    for l in amrFile.readlines():
        if "amrfinder DB:" in l.strip():
            amrfinderDB_version = l.strip().split(':')[1].strip()

with open(sys.argv[2], 'r') as krakenFile:
    for l in krakenFile.readlines():
        if "kraken DB:" in l.strip():
            krakenDB_version = l.strip().split(':')[1].strip()

files = glob.glob('*.tsv')

dfs = []

for file in files:
    df = pd.read_csv(file, header=0, delimiter='\t')
    dfs.append(df)

merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='left'), dfs)
merged = merged.assign(krakenDB=krakenDB_version)
merged = merged.assign(amrDB=amrfinderDB_version)
merged = merged.assign(spriggan=sys.argv[3])
merged = merged[['Sample','Total Reads','Reads Removed','Median Coverage','Average Coverage','Contigs','Assembly Length (bp)','N50','Primary Species (%)','Secondary Species (%)','Unclassified Reads (%)','krakenDB','MLST Scheme','Gene','Coverage','Identity','Selected AMR Genes','Selected AMR Genes Coverage','Selected AMR Genes Identity','amrDB','spriggan']]
merged = merged.rename(columns={'Contigs':'Contigs (#)','Average Coverage':'Mean Coverage','Gene':'AMR','Coverage':'AMR Coverage','Identity':'AMR Identity','krakenDB':'Kraken Database Verion','amrDB':'AMRFinderPlus Database Version','spriggan':'Spriggan Version'})

merged.to_csv('spriggan_report.csv', index=False, sep=',', encoding='utf-8')