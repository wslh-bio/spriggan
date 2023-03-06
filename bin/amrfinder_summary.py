#!/usr/bin/env python3
import os
import glob
import pandas as pd

# get amrfinder output files and set up lists
files = glob.glob('data/*.amr.tsv')
dfs = []
all_ar_dfs = []
selected_ar_dfs = []

# function for cleanining up amrfinder output
def pretty_df(data,sample):
    data.columns = data.columns.str.replace(' ', '_')
    data = data.assign(Sample=sample)
    data = data[['Sample','Gene_symbol','%_Coverage_of_reference_sequence','%_Identity_to_reference_sequence']]
    pretty_data = data.set_axis(['Sample', 'Gene', 'Coverage', 'Identity'], axis=1, inplace=False)
    return pretty_data

# function for joining amrfinder results by a delimiter
def join_df(data,sample,delim):
    gene = data['Gene'].tolist()
    gene = delim.join(gene)
    coverage = data['Coverage'].tolist()
    coverage = delim.join(map(str, coverage))
    identity = data['Identity'].tolist()
    identity = delim.join(map(str, identity))
    joined_data = [[sample,gene,coverage,identity]]
    joined_data = pd.DataFrame(joined_data, columns = ['Sample', 'Gene', 'Coverage', 'Identity'])
    return joined_data

for file in files:
    # get sample id from file name
    sample_id = os.path.basename(file).split('.')[0]
    # read in amrfinder results as data frame
    df = pd.read_csv(file, header=0, delimiter='\t')

    # clean up data frame
    df = pretty_df(df,sample_id)
    dfs.append(df)

    # summarize all results
    all_ar_df = join_df(df,sample_id,';')#!/usr/bin/env python3
import os
import glob
import pandas as pd

# get amrfinder output files and set up lists
files = glob.glob('data/*.amr.tsv')
dfs = []
all_ar_dfs = []
selected_ar_dfs = []

# function for cleanining up amrfinder output
def pretty_df(data,sample):
    data.columns = data.columns.str.replace(' ', '_')
    data = data.assign(Sample=sample)
    data = data[['Sample','Gene_symbol','%_Coverage_of_reference_sequence','%_Identity_to_reference_sequence']]
    pretty_data = data.set_axis(['Sample', 'Gene', 'Coverage', 'Identity'], axis=1, inplace=False)
    return pretty_data

# function for joining amrfinder results by a delimiter
def join_df(data,sample,delim):
    gene = data['Gene'].tolist()
    gene = delim.join(gene)
    coverage = data['Coverage'].tolist()
    coverage = delim.join(map(str, coverage))
    identity = data['Identity'].tolist()
    identity = delim.join(map(str, identity))
    joined_data = [[sample,gene,coverage,identity]]
    joined_data = pd.DataFrame(joined_data, columns = ['Sample', 'Gene', 'Coverage', 'Identity'])
    return joined_data

for file in files:
    # get sample id from file name
    sample_id = os.path.basename(file).split('.')[0]
    # read in amrfinder results as data frame
    df = pd.read_csv(file, header=0, delimiter='\t')

    # clean up data frame
    if masked_df.empty:
        masked_df = masked_df.append({'Sample' : sample_id, 'Gene' : 'None', 'Coverage' : 'None','Identity' : 'None'}, ignore_index = True)
    selected_ar_df = join_df(masked_df,sample_id,';')
    selected_ar_df = selected_ar_df.set_axis(['Sample', 'Selected AMR Genes', 'Selected AMR Genes Coverage', 'Selected AMR Genes Identity'], axis=1, inplace=False)
    selected_ar_dfs.append(selected_ar_df)

# concatenate results and write to tsv
dfs = list(dfs)
if len(dfs) > 1:
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(f'amrfinder_predictions.tsv',sep='\t', index=False, header=True, na_rep='NaN')
else:
    dfs = dfs[0]
    dfs.to_csv(f'amrfinder_predictions.tsv',sep='\t', index=False, header=True, na_rep='NaN')

# concatenate joined restults and write to tsv
all_ar_dfs = list(all_ar_dfs)
if len(dfs) > 1:
    concat_all_ar_dfs = pd.concat(all_ar_dfs)
    concat_all_ar_dfs.to_csv('amrfinder_summary.tsv',sep='\t', index=False, header=True, na_rep='NaN')
else:
    concat_all_ar_dfs = all_ar_dfs[0]
    concat_all_ar_dfs.to_csv('amrfinder_summary.tsv',sep='\t', index=False, header=True, na_rep='NaN')

# concatenate selected genes and write to tsv
selected_ar_dfs = list(selected_ar_dfs)
if len(dfs) > 1:
    concat_selected_ar_dfs = pd.concat(selected_ar_dfs)
    concat_selected_ar_dfs.to_csv('selected_ar_genes.tsv',sep='\t', index=False, header=True, na_rep='NaN')
else:
    selected_ar_dfs = selected_ar_dfs[0]
    selected_ar_dfs.to_csv('selected_ar_genes.tsv',sep='\t', index=False, header=True, na_rep='NaN')