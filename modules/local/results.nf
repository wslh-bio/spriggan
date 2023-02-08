process RESULTS {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/wslh-bioinformatics/spriggan-pandas:1.3.2"

    input:
    path("bbduk_results.tsv")
    path("coverage_stats.tsv")
    path("quast_results.tsv")
    path("mlst_results.tsv")
    path("kraken_results.tsv")
    path("amrfinder_summary.tsv")
    path("selected_ar_genes.tsv")
    path(kraken_version, stageAs:"kraken_version.yml")
    path(amrfinder_version stageAs:"amrfinder_version.yml")

    output:
    path('spriggan_report.csv')

    script:
    """
    #!/usr/bin/env python3

    import os
    import glob
    import pandas as pd
    from functools import reduce

    with open('${amrfinder_version}', 'r') as amrFile:
        for l in amrFile.readlines():
            if "amrfinder DB:" in l.strip():
                amrfinderDB_version = l.strip().split(':')[1].strip()

    with open('${kraken_version}', 'r') as krakenFile:
        for l in krakenFile.readlines():
            if "kraken DB:" in l.strip():
                krakenDB_version = l.strip().split(':')[1].strip()

    files = glob.glob('*.tsv')

    dfs = []

    for file in files:
        df = pd.read_csv(file, header=0, delimiter='\\t')
        dfs.append(df)

    merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='left'), dfs)
    merged = merged.assign(krakenDB=krakenDB_version)
    merged = merged.assign(amrDB=amrfinderDB_version)
    merged = merged[['Sample','Total Reads','Reads Removed','Median Coverage','Average Coverage','Contigs','Assembly Length (bp)','N50','Primary Species (%)','Secondary Species (%)','Unclassified Reads (%)','krakenDB','MLST Scheme','Gene','Coverage','Identity','Selected AMR Genes','Selected AMR Genes Coverage','Selected AMR Genes Identity','amrDB']]
    merged = merged.rename(columns={'Contigs':'Contigs (#)','Average Coverage':'Mean Coverage','Gene':'AMR','Coverage':'AMR Coverage','Identity':'AMR Identity','krakenDB':'Kraken Database Verion','amrDB':'AMRFinderPlus Database Version'})

    merged.to_csv('spriggan_report.csv', index=False, sep=',', encoding='utf-8')
    """
}
