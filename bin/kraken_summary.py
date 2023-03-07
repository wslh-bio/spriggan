#!/usr/bin/env python3
import os
import glob
import pandas as pd
from pandas import DataFrame

# function for summarizing kraken2 report files
def summarize_kraken(file):
    # get sample id from file name
    sample_id = os.path.basename(file).split('.')[0].replace('.kraken2.txt','')
    data = []
    # read kraken2 report file
    with open(file,'r') as inFile:
        for line in inFile:
            line = line.strip()
            sline = line.split('\t')
            # get unclassified reads result (denoted by 'unclassified') and append to data
            if sline[5] == 'unclassified':
                data.append(sline)
            # get species results (denoted by 'S') and append to data
            if sline[3] == 'S':
                data.append(sline)
    # convert data list to data frame
    data_df = DataFrame(data, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
    # remove left leading spaces from the Name column
    data_df['Name'] = data_df['Name'].str.lstrip()
    # sort data frame by percentages (largest to smallest)
    data_df['Percentage'] = pd.to_numeric(data_df['Percentage'], downcast='float')
    data_df = data_df.sort_values(by=['Percentage'], ascending=False)
    # make new data frame for unclassified reads only
    unclass = data_df[data_df['Name']=='unclassified']
    # exception for if no unclassified reads found
    if unclass.empty:
        # import pandas as pd
        lst = [['0','NA','NA','NA','NA','NA']]
        unclass = pd.DataFrame(lst, columns =['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
    # subset data frame by species
    species_df = data_df[data_df['Name']!='unclassified']
    # get first two species matches (first two largest percentages) in data frame
    species_df = species_df.head(2)
    # check if species data frame has two rows
    if len(species_df) == 0:
        # add two empty rows to species data frame
        species_df = species_df.append(pd.Series(), ignore_index=True)
        species_df = species_df.append(pd.Series(), ignore_index=True)
    if len(species_df) == 1:
        # add one empty row to species data frame
        species_df = species_df.append(pd.Series(), ignore_index=True)
    # concatenate unclassified data frame and species data frame
    df_concat = pd.concat([unclass,species_df])
    # add sample name column to concatenated data frame
    df_concat = df_concat.assign(Sample=sample_id)
    # keep only Sample Percentage and Name columns in concatenated data frame
    df_concat = df_concat[['Sample','Percentage','Name']]
    # reset index of concatenated data frame using drop parameter to avoid old index added as column
    df_concat = df_concat.reset_index(drop=True)
    # add percentage sign to unclassified column
    unclassified = str(df_concat.iloc[0]['Percentage']) + '%'
    # convert to lists
    # if primary species is nan, replace with NA
    if str(df_concat.iloc[1]['Name']) == 'nan':
        primary_species = 'NA'
    # otherwise convert to (#%)
    else:
        primary_species = df_concat.iloc[1]['Name'] + ' (' + str(df_concat.iloc[1]['Percentage']) + '%)'
    # repeat for secondary species
    if str(df_concat.iloc[2]['Name']) == 'nan':
        secondary_species = 'NA'
    else:
        secondary_species = df_concat.iloc[2]['Name'] + ' (' + str(df_concat.iloc[2]['Percentage']) + '%)'
    # list of lists
    combined = [[sample_id, unclassified, primary_species, secondary_species]]
    # convert list of lists to data frame
    combined_df = DataFrame(combined, columns=['Sample','Unclassified Reads (%)','Primary Species (%)','Secondary Species (%)'])
    return combined_df
# get all kraken2 report files
files = glob.glob("data/*.kraken2.txt")

# summarize kraken2 report files
results = map(summarize_kraken, files)

# concatenate summary results and write to tsv
results = list(results)

if len(results) > 1:
    data_concat = pd.concat(results)
    data_concat.to_csv(f'kraken_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
else:
    results = results[0]
    results.to_csv(f'kraken_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')