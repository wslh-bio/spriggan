#!/usr/bin/env python3

import glob
import pandas as pd
from pandas import DataFrame

files = glob.glob('data/*.mlst.tsv')
dfs = []
for file in files:
    df = pd.read_csv(file, sep='\\t')
    dfs.append(df)

# concatenate dfs and write data frame to file
if len(dfs) > 1:
    dfs_concat = pd.concat(dfs)
    dfs_concat['MLST Scheme'] = dfs_concat['MLST Scheme'].str.replace('-:NA', 'No Scheme Available')
    dfs_concat.to_csv(f'mlst_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
else:
    dfs = dfs[0]
    dfs['MLST Scheme'] = dfs['MLST Scheme'].str.replace('-:NA', 'No Scheme Available')
    dfs.to_csv(f'mlst_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')