#!/usr/bin/env nextflow

//Description:
//Author: Kelsey Florek and Abigail Shockey
//email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu

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

Channel
  .fromPath("$baseDir/multiqc_config.yaml")
  .set { multiqc_config }

//Preprocess reads - change names
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(outfiles) into read_files_fastqc, read_files_trimming

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

//Trim reads and remove PhiX contamination
process clean_reads {
  tag "$name"
  publishDir "${params.outdir}/trimming", mode: 'copy',pattern:"*.trim.txt"

  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_shovill, cleaned_reads_fastqc, read_files_kraken
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats
  file("${name}.trim.txt") into bbduk_files
  tuple file("${name}.phix.stats.txt"),file("${name}.adapters.stats.txt"),file("${name}.trim.txt") into multiqc_clean_reads

  script:
  """
  bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.trimmed_1.fastq.gz out2=${name}.trimmed_2.fastq.gz qtrim=${params.trimdirection} qtrim=${params.qualitytrimscore} minlength=${params.minlength} tbo tbe &> ${name}.out

  repair.sh in1=${name}.trimmed_1.fastq.gz in2=${name}.trimmed_2.fastq.gz out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz

  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo

  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_clean_1.fastq.gz out2=${name}_clean_2.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt

  grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${name}.out > ${name}.trim.txt
  """
}

//Combine raw reads channel and cleaned reads channel
combined_reads = read_files_fastqc.concat(cleaned_reads_fastqc)

//FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from combined_reads

  output:
  file("*_fastqc.{zip,html}") into fastqc_results, fastqc_multiqc

  script:
  """
  fastqc -q  ${reads}
  """
}

process fastqc_summary {
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  file(fastqc) from fastqc_results.collect()

  output:
  file("fastqc_summary.tsv") into fastqc_summary

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

//Assemble trimmed reads with Shovill and map reads back to assembly
process shovill {
  errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/assembled", mode: 'copy',pattern:"*.fa"
  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.sam"

  input:
  set val(name), file(reads) from cleaned_reads_shovill

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_ar, assembled_genomes_mlst
  tuple name, file("${name}.sam") into sam_files

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory} --outdir ./output --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv ./output/contigs.fa ${name}.contigs.fa
  bwa index ${name}.contigs.fa
  bwa mem ${name}.contigs.fa ${reads[0]} ${reads[1]} > ${name}.sam
  """
}

//Index and sort bam file then calculate coverage
process samtools {
  tag "$name"

  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.sorted.*"
  publishDir "${params.outdir}/alignments", mode: 'copy', pattern:"*.stats.txt*"
  publishDir "${params.outdir}/coverage", mode: 'copy', pattern:"*.depth.tsv*"

  input:
  set val(name), file(sam) from sam_files

  output:
  file("${name}.depth.tsv") into cov_files
  file("${name}.stats.txt") into stats_multiqc
  file("*.sorted.*")

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
  publishDir "${params.outdir}/coverage", mode: 'copy'

  input:
  file(cov) from cov_files.collect()

  output:
  file('coverage_stats.tsv') into coverage_tsv

  script:
  """
  #!/usr/bin/env python3
  import glob
  import os
  from numpy import median
  from numpy import average

  results = []

  files = glob.glob("*.depth.tsv*")
  for file in files:
    nums = []
    sid = os.path.basename(file).split('.')[0]
    with open(file,'r') as inFile:
      for line in inFile:
        nums.append(int(line.strip().split()[2]))
      med = int(median(nums))
      avg = int(average(nums))
      results.append(f"{sid}\\t{med}\\t{avg}\\n")

  with open('coverage_stats.tsv', 'w') as outFile:
    outFile.write("Sample\\tMedian Coverage\\tAverage Coverage\\n")
    for result in results:
      outFile.write(result)
  """
}

//Run Quast on assemblies
process quast {
  tag "$name"

  errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy',pattern: "${name}.quast.tsv"

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_files
  file("${name}.report.quast.tsv") into quast_multiqc

  script:
  """
  quast.py ${assembly} -o .
  mv report.tsv ${name}.report.quast.tsv
  mv transposed_report.tsv ${name}.quast.tsv
  """
}

process quast_summary {
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  file(files) from quast_files.collect()

  output:
  file("quast_results.tsv") into quast_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  files = glob.glob("*.quast.tsv")

  dfs = []

  for file in files:
      sample_id = os.path.basename(file).split(".")[0]
      df = pd.read_csv(file, sep='\\t')
      df = df.iloc[:,[1,7,17]]
      df = df.assign(Sample=sample_id)
      df = df.rename(columns={'# contigs (>= 0 bp)':'Contigs','Total length (>= 0 bp)':'Assembly Length (bp)'})
      df = df[['Sample', 'Contigs','Assembly Length (bp)', 'N50']]
      dfs.append(df)

  dfs_concat = pd.concat(dfs)
  dfs_concat.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//MLST Step 1: Run mlst
process mlst {
  tag "$name"

  publishDir "${params.outdir}/mlst", mode: 'copy', pattern: "*.mlst.tsv*"

  input:
  set val(name), file(input) from assembled_genomes_mlst

  output:
  file("*.tsv") into mlst_results

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
  """
}

//MLST Step 2: Summarize mlst results
process mlst_formatting {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/mlst",mode:'copy'

  input:
  file(mlst) from mlst_results.collect()

  output:
  file("mlst_results.tsv") into mlst_tsv

  script:
  """
  #!/usr/bin/env python3

  import glob
  import pandas as pd
  from pandas import DataFrame

  files = glob.glob('*.mlst.tsv')
  dfs = []
  for file in files:
      df = pd.read_csv(file, sep='\\t')
      dfs.append(df)
  dfs_concat = pd.concat(dfs)
  dfs_concat['MLST Scheme'] = dfs_concat['MLST Scheme'].str.replace('-:NA', 'No Scheme Available')
  dfs_concat.to_csv(f'mlst_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//Kraken Step 1: Run Kraken
process kraken {
  tag "$name"
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*.kraken2.txt*"

  input:
  set val(name), file(reads) from read_files_kraken

  output:
  tuple name, file("${name}.kraken2.txt") into kraken_files, kraken_multiqc
  file("Kraken2_DB.txt") into kraken_version

  script:
  """
  kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${name}.kraken2.txt --paired ${reads[0]} ${reads[1]}

  ls /kraken2-db/ > Kraken2_DB.txt
  """
}

//Kraken Step 2: Summarize kraken results
process kraken_summary {
  tag "$name"
  publishDir "${params.outdir}/kraken",mode:'copy'

  input:
  file(files) from kraken_files.collect()

  output:
  file("kraken_results.tsv") into kraken_tsv
  file("kraken_results.tsv") into kraken_amr

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

  # get all kraken2 report files and make results list
  files = glob.glob("*.kraken2.txt*")

  # summarize kraken2 report files
  results = map(summarize_kraken, files)

  # concatenate summary results and write to tsv
  data_concat = pd.concat(results)
  data_concat.to_csv(f'kraken_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//AR Setup amrfinder files
process amrfinder_setup {
  tag "$name"

  input:
  file(kraken) from kraken_amr
  set val(name), file(input) from assembled_genomes_ar

  output:
  tuple name, file("${name}.*.fa") into ar_input

  script:
  """
  #!/usr/bin/env python3

  import pandas as pd
  import glob
  import shutil

  species = ['Acinetobacter_baumannii','Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_aureus','Staphylococcus_pseudintermedius','Streptococcus_agalactiae','Streptococcus_pneumoniae','Streptococcus_pyogenes']
  genus = ['Campylobacter','Escherichia','Klebsiella','Salmonella']

  genomeFile = '${input}'
  sid = genomeFile.split('.')[0]

  df = pd.read_csv('kraken_results.tsv', header=0, delimiter='\\t')
  df = df[df['Sample'] == sid]
  taxa = df.iloc[0]['Primary Species (%)']
  taxa = taxa.split(' ')
  taxa_species = taxa[0] + '_' + taxa[1]
  taxa_genus = taxa[0]

  if any(x in taxa_species for x in species):
      shutil.copyfile(genomeFile, f'{sid}.{taxa_species}.fa')
  elif any(x in taxa_genus for x in genus):
      shutil.copyfile(genomeFile, f'{sid}.{taxa_genus}.fa')
  elif taxa_genus == 'Shigella':
      shutil.copyfile(genomeFile, f'{sid}.Escherichia.fa')
  else:
      shutil.copyfile(genomeFile, f'{sid}.NA.fa')
  """
}

//AR Step 2: Run amrfinder
process amrfinder {
  tag "$name"
  publishDir "${params.outdir}/amrfinder",mode: 'copy',pattern:"*.amr.tsv"

  input:
  set val(name), file(input) from ar_input

  output:
  tuple name, file("${name}.amr.tsv") into ar_predictions
  file("AMRFinderPlus_DB.txt") into amrfinder_version

  script:
  """
  #!/usr/bin/env python3

  import subprocess as sub
  import shlex
  import glob
  import shutil

  organisms = ['Acinetobacter_baumannii','Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_aureus','Staphylococcus_pseudintermedius','Streptococcus_agalactiae','Streptococcus_pneumoniae','Streptococcus_pyogenes','Campylobacter','Escherichia','Klebsiella','Salmonella','Escherichia']

  fastaFile = '${input}'
  sid = fastaFile.split('.')[0]
  organism = fastaFile.split('.')[1]

  if any(x in organism for x in organisms):
      outFile = open(f'{sid}.amr.tsv','w')
      cmd = shlex.split(f'amrfinder -n {sid}.{organism}.fa --organism {organism}')
      sub.Popen(cmd, stdout=outFile).wait()
  else:
      outFile = open(f'{sid}.amr.tsv','w')
      cmd = shlex.split(f'amrfinder -n {sid}.{organism}.fa')
      sub.Popen(cmd, stdout=outFile).wait()

  versionFile = "/amrfinder/data/latest/version.txt"
  shutil.copy(versionFile,"AMRFinderPlus_DB.txt")
  """
}

//AR Step 3: Summarize amrfinder+ results
process amrfinder_summary {
  publishDir "${params.outdir}/amrfinder",mode:'copy'

  input:
  file(predictions) from ar_predictions.collect()

  output:
  file("ar_predictions.tsv")
  file("ar_summary.tsv") into ar_tsv
  file("important_ar_genes.tsv") into important_ar_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd

  files = glob.glob('*.amr.tsv')
  dfs = []
  all_ar_dfs = []
  important_ar_dfs = []

  def pretty_df(data,sample):
      data.columns = data.columns.str.replace(' ', '_')
      data = data.assign(Sample=sample)
      data = data[['Sample','Gene_symbol','%_Coverage_of_reference_sequence','%_Identity_to_reference_sequence']]
      pretty_data = data.set_axis(['Sample', 'Gene', 'Coverage', 'Identity'], axis=1, inplace=False)
      return pretty_data

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
      sample_id = os.path.basename(file).split('.')[0]
      df = pd.read_csv(file, header=0, delimiter='\\t')

      df = pretty_df(df,sample_id)
      dfs.append(df)

      all_ar_df = join_df(df,sample_id,';')
      all_ar_dfs.append(all_ar_df)

      mask = df['Gene'].str.contains(${params.important_genes}, case=False, na=False)
      masked_df = df[mask]
      important_ar_df = join_df(masked_df,sample_id,' ')
      important_ar_df = important_ar_df.set_axis(['Sample', 'Important Genes', 'Important Genes Coverage', 'Important Genes Identity'], axis=1, inplace=False)
      important_ar_dfs.append(important_ar_df)

  concat_dfs = pd.concat(dfs)
  concat_dfs.to_csv('ar_predictions.tsv',sep='\\t', index=False, header=True, na_rep='NaN')

  concat_all_ar_dfs = pd.concat(all_ar_dfs)
  concat_important_ar_dfs = pd.concat(important_ar_dfs)

  concat_all_ar_dfs.to_csv('ar_summary.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  concat_important_ar_dfs.to_csv('important_ar_genes.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

process bbduk_summary {
  publishDir "${params.outdir}/trimming",mode:'copy'

  input:
  file(files) from bbduk_files.collect()

  output:
  file("bbduk_results.tsv") into bbduk_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  files = glob.glob("*.txt")

  results = []
  for file in files:
      sample_id = os.path.basename(file).split(".")[0]
      vals = []
      vals.append(sample_id)
      with open(file,"r") as inFile:
          for i, line in enumerate(inFile):
              if i == 0:
                  num_reads = line.strip().split("\\t")[1].replace(" reads ","")
                  vals.append(num_reads)
              if i == 3:
                  rm_reads = line.strip().split("\\t")[1].replace("reads ","")
                  rm_reads = rm_reads.rstrip()
                  vals.append(rm_reads)
      results.append(vals)

  df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])
  df.to_csv(f'bbduk_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//Merge results
process merge_results {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  file(bbduk) from bbduk_tsv
  file(quast) from quast_tsv
  file(coverage) from coverage_tsv
  file(mlst) from mlst_tsv
  file(kraken) from kraken_tsv
  file(amr) from ar_tsv
  file(important_ar) from important_ar_tsv
  file(vkraken) from kraken_version.first()
  file(vamrfinder) from amrfinder_version.first()

  output:
  file('spriggan_report.csv')

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
  merged = merged[['Sample','Total Reads','Reads Removed','Median Coverage','Average Coverage','Contigs','Assembly Length (bp)','N50','Primary Species (%)','Secondary Species (%)','Unclassified Reads (%)','krakenDB','MLST Scheme','Important Genes','Gene','Coverage','Identity','Important Genes Coverage','Important Genes Identity','amrDB']]
  merged = merged.rename(columns={'Contigs':'Contigs (#)','Average Coverage':'Mean Coverage','Gene':'AMR','Coverage':'AMR Coverage','Identity':'AMR Identity','krakenDB':'Kraken Database Verion','amrDB':'AMRFinderPlus Database Version'})

  merged.to_csv('spriggan_report.csv', index=False, sep=',', encoding='utf-8')
  """
}

//QC Step: MultiQC
process multiqc {
  publishDir "${params.outdir}",mode:'copy'

  input:
  file(a) from multiqc_clean_reads.collect()
  file(b) from fastqc_multiqc.collect()
  file(c) from stats_multiqc.collect()
  file(d) from kraken_multiqc.collect()
  file(e) from quast_multiqc.collect()
  file(config) from multiqc_config

  output:
  file("*.html") into multiqc_output

  script:
  """
  multiqc -c ${config} .
  """
}
