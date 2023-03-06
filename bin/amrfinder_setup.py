#!/usr/bin/env python3
import pandas as pd
import shutil

# species and genus lists
species = ['Acinetobacter_baumannii','Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_aureus','Staphylococcus_pseudintermedius','Streptococcus_agalactiae','Streptococcus_pneumoniae','Streptococcus_pyogenes']
genus = ['Escherichia','Salmonella']

# get sample name from fasta file
genomeFile = '${assembly}'
sid = genomeFile.split('.')[0]

# read in kraken results as data frame
df = pd.read_csv('kraken_results.tsv', header=0, delimiter='\\t')

# subset data frame by sample id
df = df[df['Sample'] == sid]

# get primary species and genus identified
if df.empty:
    taxa_species = 'NA'
    taxa_genus = 'NA'
else:
    taxa = df.iloc[0]['Primary Species (%)']
    taxa = taxa.split(' ')
    taxa_species = taxa[0] + '_' + taxa[1]
    taxa_genus = taxa[0]

# add taxa or genus name to file name if present in lists
if any(x in taxa_species for x in species):
    shutil.copyfile(genomeFile, f'{sid}.{taxa_species}.fa')
elif taxa_species == 'Campylobacter_coli' or taxa_species == 'Campylobacter_jejuni':
    shutil.copyfile(genomeFile, f'{sid}.Campylobacter.fa')
elif taxa_species == 'Enterococcus_hirae':
    shutil.copyfile(genomeFile, f'{sid}.Enterococcus_faecium.fa')
elif taxa_genus == 'Shigella':
    shutil.copyfile(genomeFile, f'{sid}.Escherichia.fa')
elif taxa_species == 'Klebsiella_aerogenes':
    shutil.copyfile(genomeFile, f'{sid}.Klebsiella.fa')
elif taxa_species == 'Neisseria_gonorrhoeae' or taxa_species == 'Neisseria_meningitidis':
    shutil.copyfile(genomeFile, f'{sid}.Neisseria.fa')
elif taxa_species == 'Streptococcus_mitis':
        shutil.copyfile(genomeFile, f'{sid}.Streptococcus_pneumoniae.fa')
elif any(x in taxa_genus for x in genus):
    shutil.copyfile(genomeFile, f'{sid}.{taxa_genus}.fa')
else:
    shutil.copyfile(genomeFile, f'{sid}.NA.fa')