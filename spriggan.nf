#!/usr/bin/env nextflow

//Description:
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

nextflow.enable.dsl=2

params.test = false

if(params.test){
  testIDS = ['SRR14131356','SRR14131352','SRR14311556','SRR14568713',
  'SRR14131354','SRR14613517','SRR14613503','SRR14613708','SRR14613700',
  'SRR14616016']

  println "Running test analysis using the following samples:"
  println testIDS
  Channel
      .fromSRA(testIDS)
      .set { raw_reads }

} else{
  //setup channel to read in and pair the fastq files
  Channel
      .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
      .set { raw_reads }
}

//Preprocessing Step: Change read names
process preProcess {
  //publishDir "${params.outdir}/reads", mode: 'copy', pattern:"*.gz"

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("*_R{1,2}.fastq.gz"), emit: processed_reads

  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    outfiles = ["${name}_R1.fastq.gz","${name}_R2.fastq.gz"]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
    outfiles = reads
    """
    """
  }
}

//QC Step: Trim reads and remove adapters and remove PhiX contamination
process clean_reads {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/trimming/stats", mode: 'copy', pattern:"*.trim.txt"
  publishDir "${params.outdir}/trimming/reads", mode: 'copy', pattern:"*.gz"

  input:
  tuple val(name), path(processed_reads)

  output:
  tuple val(name), path("${name}_clean{_1,_2}.fastq.gz"), emit: cleaned_reads
  path("${name}.trim.txt"), emit: bbduk_files
  path("${name}.adapter.stats.txt"), emit: multiqc_adapters

  script:
  """
  bbduk.sh in1=${processed_reads[0]} in2=${processed_reads[1]} out1=${name}_clean_1.fastq.gz out2=${name}_clean_2.fastq.gz outm=${name}.adapters.fq qtrim=${params.trimdirection} trimq=${params.qualitytrimscore} minlength=${params.minlength} ref=/bbmap/resources/adapters.fa stats=${name}.adapter.stats.txt k=31 hdist=1 tpe tbo &> ${name}.out
  grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${name}.out > ${name}.trim.txt
  """
}

//Summary Step: Summarize BBDuk results
process bbduk_summary {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/trimming",mode:'copy'

  input:
  path("data*/*")

  output:
  path("bbduk_results.tsv"), emit: bbduk_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  # function for summarizing bbduk output
  def summarize_bbduk(file):
      # get sample id from file name and set up data list
      sample_id = os.path.basename(file).split(".")[0]
      data = []
      data.append(sample_id)
      with open(file,"r") as inFile:
          for i, line in enumerate(inFile):
              # get total number of reads
              if i == 0:
                  num_reads = line.strip().split("\\t")[1].replace(" reads ","")
                  data.append(num_reads)
              # get total number of reads removed
              if i == 3:
                  rm_reads = line.strip().split("\\t")[1].replace("reads ","")
                  rm_reads = rm_reads.rstrip()
                  data.append(rm_reads)
      return data

  # get all bbduk output files
  files = glob.glob("data*/*.trim.txt")

  # summarize bbduk output files
  results = map(summarize_bbduk,files)

  # convert results to data frame and write to tsv
  df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])
  df.to_csv(f'bbduk_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}


//QC Step: Run FastQC on processed and cleaned reads
process fastqc {
  //errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  tuple val(name), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
}

//Summary Step: Summarize FastQC results
process fastqc_summary {
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  path(fastqc_results)

  output:
  path("fastqc_summary.tsv"), emit: fastqc_summary

  shell:
  """
  zips=`ls *.zip`

  for i in \$zips; do
      unzip -o \$i &>/dev/null;
  done

  fq_folders=\${zips}

  for folder in \$fq_folders; do
    folder=\${folder%.*}
    cat \$folder/summary.txt >> fastqc_summary.tsv
    ls .
  done;

  sed -i 's/.fastq.gz//g' fastqc_summary.tsv
  """
}

//Assembly step: Assemble trimmed reads with Shovill and map reads back to assembly
process shovill {
  //errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/assembled", mode: 'copy',pattern:"*.fa"
  publishDir "${params.outdir}/mapping/sams", mode: 'copy', pattern:"*.sam"

  input:
  tuple val(name), path(cleaned_reads)

  output:
  tuple val(name), path("${name}.contigs.fa"), emit: assembled_genomes
  tuple val(name), path("${name}.sam"), emit: sam_files

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory} --outdir ./output --R1 ${cleaned_reads[0]} --R2 ${cleaned_reads[1]} --force
  mv ./output/contigs.fa ${name}.contigs.fa
  bwa index ${name}.contigs.fa
  bwa mem ${name}.contigs.fa ${cleaned_reads[0]} ${cleaned_reads[1]} > ${name}.sam
  """
}

//Index and sort bam file then calculate coverage
process samtools {
  //errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/mapping/bams", mode: 'copy', pattern:"*.sorted.bam*"
  publishDir "${params.outdir}/mapping/depth", mode: 'copy', pattern:"*.depth.tsv"
  publishDir "${params.outdir}/mapping/stats", mode: 'copy', pattern:"*.stats.txt"

  input:
  tuple val(name), path(sam_files)

  output:
  path("*.depth.tsv"), emit: cov_files
  path("*.stats.txt"), emit: stats_multiqc
  path("*.sorted.*")

  shell:
  """
  samtools view -S -b ${name}.sam > ${name}.bam
  samtools sort ${name}.bam > ${name}.sorted.bam
  samtools index ${name}.sorted.bam
  samtools depth -a ${name}.sorted.bam > ${name}.depth.tsv
  samtools stats ${name}.sorted.bam > ${name}.stats.txt
  """
}

//Calculate coverage stats
process coverage_stats {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/mapping", mode: 'copy', pattern:"*.tsv"

  input:
  path("data*/*")

  output:
  path('coverage_stats.tsv'), emit: coverage_tsv

  script:
  """
  #!/usr/bin/env python3
  import glob
  import os
  from numpy import median
  from numpy import average

  # function for summarizing samtools depth files
  def summarize_depth(file):
      # get sample id from file name and set up data list
      sid = os.path.basename(file).split('.')[0]
      data = []
      # open samtools depth file and get depth
      with open(file,'r') as inFile:
          for line in inFile:
              data.append(int(line.strip().split()[2]))
      # get median and average depth
      med = int(median(data))
      avg = int(average(data))
      # return sample id, median and average depth
      result = f"{sid}\\t{med}\\t{avg}\\n"
      return result

  # get all samtools depth files
  files = glob.glob("data*/*.depth.tsv")
  print(files)

  # summarize samtools depth files
  results = map(summarize_depth,files)

  # write results to file
  with open('coverage_stats.tsv', 'w') as outFile:
      outFile.write("Sample\\tMedian Coverage\\tAverage Coverage\\n")
      for result in results:
          outFile.write(result)
  """
}

//QC Step: Run QUAST on assemblies
process quast {
//  errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/quast",mode:'copy',pattern: "*.quast.tsv"

  input:
  tuple val(name), path(assembled_genomes)

  output:
  path("*.transposed.quast.tsv"), emit: quast_files
  path("*.report.quast.tsv"), emit: multiqc_quast

  script:
  """
  quast.py ${name}.contigs.fa -o .
  mv report.tsv ${name}.report.quast.tsv
  mv transposed_report.tsv ${name}.transposed.quast.tsv
  """
}

//QC Step: Run QUAST on assemblies
process quast_summary {
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  path("data*/*")

  output:
  path("quast_results.tsv"), emit: quast_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  # function for summarizing quast output
  def summarize_quast(file):
      # get sample id from file name and set up data list
      sample_id = os.path.basename(file).split(".")[0]
      # read in data frame from file
      df = pd.read_csv(file, sep='\\t')
      # get contigs, total length and assembly length columns
      df = df.iloc[:,[1,7,17]]
      # assign sample id as column
      df = df.assign(Sample=sample_id)
      # rename columns
      df = df.rename(columns={'# contigs (>= 0 bp)':'Contigs','Total length (>= 0 bp)':'Assembly Length (bp)'})
      # re-order data frame
      df = df[['Sample', 'Contigs','Assembly Length (bp)', 'N50']]
      return df

  # get quast output files
  files = glob.glob("data*/*.transposed.quast.tsv")

  # summarize quast output files
  dfs = map(summarize_quast,files)

  # concatenate dfs and write data frame to file
  dfs_concat = pd.concat(dfs)
  dfs_concat.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//MLST Step 1: Run mlst
process mlst {
  errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/mlst/schemes", mode: 'copy', pattern: "*.mlst.tsv"
  publishDir "${params.outdir}/mlst/alleles", mode: 'copy', pattern: "*.alleles.tsv"

  input:
  tuple val(name), path(input)

  output:
  path("*.mlst.tsv"), emit: mlst_files
  path("*.alleles.tsv")

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
  """
}

//MLST Step 2: Summarize mlst results
process mlst_summary {
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/mlst",mode:'copy'

  input:
  path("data*/*")

  output:
  path("mlst_results.tsv"), emit: mlst_tsv

  script:
  """
  #!/usr/bin/env python3

  import glob
  import pandas as pd
  from pandas import DataFrame

  files = glob.glob('data*/*.mlst.tsv')
  dfs = []
  for file in files:
      df = pd.read_csv(file, sep='\\t')
      dfs.append(df)
  dfs_concat = pd.concat(dfs)
  dfs_concat['MLST Scheme'] = dfs_concat['MLST Scheme'].str.replace('-:NA', 'No Scheme Available')
  dfs_concat.to_csv(f'mlst_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//Classification Step: Run Kraken
process kraken {
  tag "$name"
  //errorStrategy 'ignore'
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*.kraken2.txt*"

  input:
  tuple val(name), path(cleaned_reads)

  output:
  path("${name}.kraken2.txt"), emit: kraken_files
  path("Kraken2_DB.txt"), emit: kraken_version

  script:
  """
  kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${name}.kraken2.txt --paired ${cleaned_reads[0]} ${cleaned_reads[1]}

  ls /kraken2-db/ > Kraken2_DB.txt
  """
}

//Summary Step: Summarize kraken results
process kraken_summary {
  tag "$name"
//   errorStrategy 'ignore'
  publishDir "${params.outdir}/kraken",mode:'copy'

  input:
  path("data*/*")

  output:
  path("kraken_results.tsv"), emit: kraken_tsv

  script:
  """
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
              sline = line.split('\\t')
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
      unclassified = df_concat.iloc[0]['Percentage'] + '%'
      # convert to lists
      # if primary species is nan, replace with NA
      if str(df_concat.iloc[1]['Name']) == 'nan':
          primary_species = 'NA'
      # otherwise convert to (#%)
      else:
          primary_species = df_concat.iloc[1]['Name'] + ' (' + df_concat.iloc[1]['Percentage'] + '%)'
      # repeat for secondary species
      if str(df_concat.iloc[2]['Name']) == 'nan':
          secondary_species = 'NA'
      else:
          secondary_species = df_concat.iloc[2]['Name'] + ' (' + df_concat.iloc[2]['Percentage'] + '%)'
      # list of lists
      combined = [[sample_id, unclassified, primary_species, secondary_species]]
      # convert list of lists to data frame
      combined_df = DataFrame(combined, columns=['Sample','Unclassified Reads (%)','Primary Species (%)','Secondary Species (%)'])
      return combined_df
  # get all kraken2 report files
  files = glob.glob("data*/*.kraken2.txt")
  # summarize kraken2 report files
  results = map(summarize_kraken, files)
  # concatenate summary results and write to tsv
  data_concat = pd.concat(results)
  data_concat.to_csv(f'kraken_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//AR Setup amrfinder files
process amrfinder_setup {
  //errorStrategy 'ignore'
  tag "$name"

  input:
  path("kraken_results.tsv")
  tuple val(name), path(assembly)

  output:
  tuple val(name), path("${name}.*.fa"), emit: amrfinder_input

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import glob
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
  """
}

//AR Step 2: Run amrfinder
process amrfinder {
  //errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/amrfinder",mode: 'copy',pattern:"*.amr.tsv"

  input:
  tuple val(name), path(amrfinder_input)

  output:
  path("${name}.amr.tsv"), emit: amrfinder_predictions
  path("AMRFinderPlus_DB.txt"), emit: amrfinder_version

  script:
  """
  #!/usr/bin/env python3

  import subprocess as sub
  import shlex
  import glob
  import shutil

  # organism list
  organisms = ['Acinetobacter_baumannii','Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_aureus','Staphylococcus_pseudintermedius','Streptococcus_agalactiae','Streptococcus_pneumoniae','Streptococcus_pyogenes','Campylobacter','Escherichia','Klebsiella','Salmonella','Escherichia']

  # get sample id and organism name from fasta file
  fastaFile = '${amrfinder_input}'
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

  # get version information from version file
  versionFile = "/amrfinder/data/latest/version.txt"
  shutil.copy(versionFile,"AMRFinderPlus_DB.txt")
  """
}

//AR Step 3: Summarize amrfinder+ results
process amrfinder_summary {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/amrfinder",mode:'copy'

  input:
  path("data*/*")

  output:
  path("amrfinder_predictions.tsv")
  path("amrfinder_summary.tsv"), emit: amrfinder_tsv
  path("selected_ar_genes.tsv"), emit: selected_ar_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd

  # get amrfinder output files and set up lists
  files = glob.glob('data*/*.amr.tsv')
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
      df = pd.read_csv(file, header=0, delimiter='\\t')

      # clean up data frame
      df = pretty_df(df,sample_id)
      dfs.append(df)

      # summarize all results
      all_ar_df = join_df(df,sample_id,';')
      all_ar_dfs.append(all_ar_df)

      # subset data frame by selected genes
      mask = df['Gene'].str.contains(${params.selected_genes}, case=False, na=False)
      masked_df = df[mask]
      # check if any select genes were found
      if masked_df.empty:
          masked_df = masked_df.append({'Sample' : sample_id, 'Gene' : 'None', 'Coverage' : 'None','Identity' : 'None'}, ignore_index = True)
      selected_ar_df = join_df(masked_df,sample_id,';')
      selected_ar_df = selected_ar_df.set_axis(['Sample', 'Selected AMR Genes', 'Selected AMR Genes Coverage', 'Selected AMR Genes Identity'], axis=1, inplace=False)
      selected_ar_dfs.append(selected_ar_df)

  # concatenate results and write to tsv
  concat_dfs = pd.concat(dfs)
  concat_dfs.to_csv('amrfinder_predictions.tsv',sep='\\t', index=False, header=True, na_rep='NaN')

  # concatenate joined restults and write to tsv
  concat_all_ar_dfs = pd.concat(all_ar_dfs)
  concat_selected_ar_dfs = pd.concat(selected_ar_dfs)

  # concatenate selected genes and write to tsv
  concat_all_ar_dfs.to_csv('amrfinder_summary.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  concat_selected_ar_dfs.to_csv('selected_ar_genes.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//Merge results
process merge_results {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  path("bbduk_results.tsv")
  path("coverage_stats.tsv")
  path("quast_results.tsv")
  path("mlst_results.tsv")
  path("kraken_results.tsv")
  path("amrfinder_summary.tsv")
  path("selected_ar_genes.tsv")
  path("Kraken2_DB.txt")
  path("AMRFinderPlus_DB.txt")

  output:
  path('spriggan_report.csv')

  script:
  """
  #!/usr/bin/env python3

  import os
  import glob
  import pandas as pd
  from functools import reduce

  with open('AMRFinderPlus_DB.txt', 'r') as amrFile:
      amrfinderDB_version = amrFile.readline().strip()

  with open('Kraken2_DB.txt', 'r') as krakenFile:
      krakenDB_version = krakenFile.readline().strip()

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

Channel
  .fromPath("$baseDir/configs/multiqc_config.yaml")
  .set { multiqc_config }

//QC Step: MultiQC
process multiqc {
  publishDir "${params.outdir}",mode:'copy'

  input:
  path("data*/*")
  path(config)

  output:
  path("*.html"), emit: multiqc_output

  script:
  """
  multiqc -c ${config} .
  """
}

workflow {

    preProcess(raw_reads)

    clean_reads(preProcess.out.processed_reads)

    bbduk_summary(clean_reads.out.bbduk_files.collect())

    processed = preProcess.out.processed_reads
    cleaned = clean_reads.out.cleaned_reads
    combined_reads = processed.concat(cleaned)

    fastqc(combined_reads)

    fastqc_summary(fastqc.out.collect())

    shovill(clean_reads.out.cleaned_reads)

    samtools(shovill.out.sam_files)

    coverage_stats(samtools.out.cov_files.collect())

    quast(shovill.out.assembled_genomes)

    quast_summary(quast.out.quast_files.collect())

    mlst(shovill.out.assembled_genomes)

    mlst_summary(mlst.out.mlst_files.collect())

    kraken(clean_reads.out.cleaned_reads)

    kraken_summary(kraken.out.kraken_files.collect())

    amrfinder_setup(kraken_summary.out.kraken_tsv,shovill.out.assembled_genomes)

    amrfinder(amrfinder_setup.out.amrfinder_input)

    amrfinder_summary(amrfinder.out.amrfinder_predictions.collect())

    merge_results(bbduk_summary.out.bbduk_tsv,coverage_stats.out.coverage_tsv,quast_summary.out.quast_tsv,mlst_summary.out.mlst_tsv,kraken_summary.out.kraken_tsv,amrfinder_summary.out.amrfinder_tsv,amrfinder_summary.out.selected_ar_tsv,kraken.out.kraken_version.first(),amrfinder.out.amrfinder_version.first())

    multiqc(clean_reads.out.bbduk_files.mix(clean_reads.out.bbduk_files,clean_reads.out.multiqc_adapters,fastqc.out.fastqc_results,samtools.out.stats_multiqc,kraken.out.kraken_files,quast.out.multiqc_quast).collect(),multiqc_config)
}
