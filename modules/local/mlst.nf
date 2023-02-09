process MLST {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/wslh-bioinformatics/spriggan-mlst:2.19.0"

    input:
    tuple val(meta), path(input)

    output:
    path("*.mlst.tsv")     , emit: mlst_files
    path "versions.yml"    , emit: versions
    path("*.alleles.tsv")


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import os
    import subprocess as sub
    import shlex
    import pandas as pd
    from functools import reduce
    import shutil
    import glob

    # Function for running mlst on samples with multiple schemes
    def run_schemes(first_scheme,all_schemes,sample):
        # Subtract selected scheme from list of schemes to get remaining schemes
        remaining_schemes = list(set(all_schemes) - set([first_scheme]))
        # Run mlst on remaining schemes
        for i in remaining_schemes:
            outFile = open(f'{sample}.{i}.tsv','w')
            cmd = shlex.split(f'mlst --nopath --exclude sthermophilus --scheme {i} {sample}.contigs.fa')
            sub.Popen(cmd, stdout=outFile).wait()

    # Read in fasta file
    inFile = '${input}'
    sid = inFile.split('.')[0]

    # Open outfile and run mlst
    outFile = open(f'{sid}.tsv','w')
    cmd = shlex.split(f'mlst --nopath --exclude sthermophilus {sid}.contigs.fa')
    sub.Popen(cmd, stdout=outFile).wait()

    # Lists organisms with multiple schemes
    abaumannii_schemes = ['abaumannii','abaumannii_2']
    ecoli_schemes = ['ecoli','ecoli_2']
    leptospira_schemes = ['leptospira','leptospira_2','leptospira_3']
    vcholerae_schemes = ['vcholerae','vcholerae2']

    # Dictionary of scheme names
    ids = {'mlstID':['abaumannii','abaumannii_2','ecoli','ecoli_2','leptospira','leptospira_2','leptospira_3','vcholerae2','vcholerae'],
    'PubMLSTID':['Oxford','Pasteur','Achtman','Pasteur ','Scheme 1','Scheme 2', 'Scheme 3','O1 and O139','']}
    ids = dict(zip(ids['mlstID'], ids['PubMLSTID']))

    # read in mlst output and get scheme
    df = pd.read_csv(f'{sid}.tsv', header=None, delimiter='\\t')
    scheme = df.iloc[0][1]

    # Add scheme to mlst file name
    if scheme == '-':
        shutil.move(f'{sid}.tsv',f'{sid}.NA.tsv')
    else:
        shutil.move(f'{sid}.tsv', f'{sid}.{scheme}.tsv')

    dfs = []

    # Check and run multiple schemes
    if any(x in scheme for x in abaumannii_schemes):
        run_schemes(scheme,abaumannii_schemes,sid)
    if any(x in scheme for x in ecoli_schemes):
        run_schemes(scheme,ecoli_schemes,sid)
    if any(x in scheme for x in leptospira_schemes):
        run_schemes(scheme,leptospira_schemes,sid)
    if any(x in scheme for x in vcholerae_schemes):
        run_schemes(scheme,vcholerae_schemes,sid)

    # Get list of mlst files and set up empty list
    mlst_files = glob.glob('*.tsv')

    # Reformat MLST results and append to empty list
    for file in mlst_files:
        df = pd.read_csv(file, header=None, delimiter='\\t')
        df[0] = df[0].str.replace('.contigs.fa', '')
        df[2] = 'ST' + df[2].astype(str)
        df[2] = df[2].str.replace('ST-', 'NA')

        if len(mlst_files) > 1:
            # Replace mlst scheme names with PubMLST scheme names
            for old, new in ids.items():
                df[1] = df[1].replace(to_replace=old, value=new)
        else:
            # Remove scheme name
            df.iloc[0,1] = ''
            df[2] = df[2].str.replace('NA', 'No scheme available')
        # Join ST to PubMLST scheme names
        df['MLST Scheme'] = df[[1,2]].agg(' '.join, axis=1)
        df = df[[0,'MLST Scheme']]
        df.columns =['Sample','MLST Scheme']
        df['MLST Scheme'] = df['MLST Scheme'].replace('\\s+', ' ', regex=True)
        df['MLST Scheme'] = df['MLST Scheme'].str.replace('NA -', 'NA', regex=True)
        dfs.append(df)

    # Merge multiple dataframes (separated by ;) and write to file
    if len(dfs) > 1:
        merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'], how='left'), dfs)
        merged = merged.reindex(sorted(merged.columns,reverse=True), axis=1)
        merged['MLST Scheme'] = merged.iloc[: , 1:].agg(';'.join, axis=1)
        merged = merged[['Sample','MLST Scheme']]
        merged.to_csv(f'{sid}.mlst.tsv', index=False, sep='\\t', encoding='utf-8')

    else:
        df = dfs[0]
        df['MLST Scheme'] = df['MLST Scheme'].str.replace(' S', 'S')
        df['MLST Scheme'] = df['MLST Scheme'].str.replace(' N', 'N')
        df.to_csv(f'{sid}.mlst.tsv', index=False, sep='\\t', encoding='utf-8')

    # rename raw mlst output (scheme files) with scheme name and add "alleles" to file name
    mlst_files = glob.glob('*.tsv')
    scheme_files = [file for file in mlst_files if not ('mlst' in file)]
    # scheme_files = [file for file in mlst_files if not ('mlst' in file or 'NA' in file)]
    for file in scheme_files:
        sid = file.split('.')[0]
        scheme = file.split('.')[1]
        shutil.move(f'{file}', f'{sid}.{scheme}.alleles.tsv')


    cmd = '''
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$(echo \$(mlst --version 2>&1) | sed 's/^.*mlst //')
    END_VERSIONS
    '''
    cmd = shlex.split(cmd)
    sub.Popen(cmd).wait()  
    """
}
