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
merged = merged[['Sample',
                 'Total Reads',
                 'Reads Removed',
                 'Median Coverage',
                 'Average Coverage',
                 'Contigs',
                 'Assembly Length (bp)',
                 'N50','Primary Species (%)',
                 'Secondary Species (%)',
                 'Unclassified Reads (%)',
                 'krakenDB','MLST Scheme',
                 'Gene','Coverage',
                 'Identity','Selected AMR Genes',
                 'Selected AMR Genes Coverage',
                 'Selected AMR Genes Identity',
                 'Expected Genome Length',
                 'Genome Length Ratio (Actual/Expected)',
                 'Species GC Content (Mean)',
                 'Sample GC Content (%)']]
merged = merged.rename(columns={'Contigs':'Contigs (#)',
                                'Average Coverage':'Mean Coverage',
                                'Gene':'AMR','Coverage':'AMR Coverage',
                                'Identity':'AMR Identity',
                                'krakenDB':'Kraken Database Verion',
                                'amrDB':'AMRFinderPlus Database Version',
                                'spriggan':'Spriggan Version'})

# Modify the MLST Scheme column to format MLST<WGS result>_<scheme used>_<organism name>

def modify_mlst_scheme(row):
    scheme_str = row['MLST Scheme']
    species = row['Primary Species (%)'].split(' (')[0]
    
    # handle missing or empty MLST Scheme
    if pd.isna(scheme_str) or scheme_str.strip() == "":
        return scheme_str
    
    # split on semicolon if present (for organisms with multiple schemes)
    schemes = scheme_str.split(';') if ';' in scheme_str else [scheme_str]
    
    # modify the MLST scheme with species appended to each
    modified = [f"{s}_{species}" for s in schemes]
    
    # join with "_" separator
    return '_'.join(modified)

merged['MLST Scheme'] = merged.apply(modify_mlst_scheme, axis=1)

merged.to_csv('spriggan_report.csv', index=False, sep=',', encoding='utf-8')