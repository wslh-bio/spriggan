#!/usr/bin/env python3

import sys
import subprocess as sub
import shlex

# organism list
organisms = ['Acinetobacter_baumannii','Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_aureus','Staphylococcus_pseudintermedius','Streptococcus_agalactiae','Streptococcus_pneumoniae','Streptococcus_pyogenes','Campylobacter','Escherichia','Klebsiella','Salmonella','Escherichia']

# get sample id and organism name from fasta file
fastaFile = sys.argv[1] 
sid = fastaFile.split('.')[0]
organism = fastaFile.split('.')[1]

# run amrfinder using --organism if present in organism list
if any(x in organism for x in organisms):
    outFile = open(f'{sid}.amr.tsv','w')
    cmd = shlex.split(f'amrfinder -n {sid}.{organism}.fa --organism {organism}')
    sub.Popen(cmd, stdout=outFile).wait()
# otherwise run amrfinder without --organism
else:
    outFile = open(f'{sid}.amr.tsv','w')
    cmd = shlex.split(f'amrfinder -n {sid}.{organism}.fa')
    sub.Popen(cmd, stdout=outFile).wait()