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
  bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.trimmed_1.fastq.gz out2=${name}.trimmed_2.fastq.gz qtrim=window,${params.windowsize} trimq=${params.qualitytrimscore} minlength=${params.minlength} tbo tbe &> ${name}.out

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
  file("fq_summary.txt") into fastqc_summary

  shell:
  """
  zips=`ls *.zip`

  for i in \$zips; do
      unzip -o \$i &>/dev/null;
  done

  fq_folders=\${zips}

  for folder in \$fq_folders; do
    folder=\${folder%.*}
    cat \$folder/summary.txt >> fq_summary.txt
    ls .
  done;

  sed -i 's/.fastq.gz//g' fq_summary.txt
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
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_annotation, assembled_genomes_ar, assembled_genomes_mash, assembled_genomes_mlst
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

  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.bam"
  publishDir "${params.outdir}/coverage", mode: 'copy', pattern:"*_depth.tsv*"

  input:
  set val(name), file(sam) from sam_files

  output:
  file("${name}_depth.tsv") into cov_files
  file("${name}.stats.txt") into stats_multiqc

  shell:
  """
  samtools view -S -b ${name}.sam > ${name}.bam
  samtools sort ${name}.bam > ${name}.sorted.bam
  samtools index ${name}.sorted.bam
  samtools depth -a ${name}.sorted.bam > ${name}_depth.tsv
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

  files = glob.glob("*_depth.tsv*")
  for file in files:
    nums = []
    sid = os.path.basename(file).split('_')[0]
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
  errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_files, quast_multiqc

  script:
  """
  quast.py ${assembly} -o .
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
  errorStrategy 'ignore'
  publishDir "${params.outdir}/mlst",mode:'copy'

  input:
  file(assemblies) from assembled_genomes_mlst.collect()

  output:
  file("mlst.tsv") into mlst_results

  script:
  """
  mlst --nopath *.fa > mlst.tsv
  """
}

//MLST Step 2: Summarize mlst results
process mlst_formatting {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/mlst",mode:'copy'

  input:
  file(mlst) from mlst_results

  output:
  file("*.tsv") into mlst_tsv

  script:
  """
  #!/usr/bin/env python3
  import csv

  string_map = {}

  with open('mlst.tsv','r') as csvfile:
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
    csvfile.seek(0)
    reader = csv.reader(csvfile,dialect,delimiter='\t')
    for row in reader:
      id_string = row[0]
      sp_string = row[1]
      st_string = row[2]
      string_map[id_string] = [sp_string,st_string]

  mlst = []
  for key in string_map:
    id = key
    id = id.replace('.contigs.fa','')
    species = string_map[key][0]
    st = string_map[key][1]
    if species == 'abaumannii':
        st = 'PubMLST ST' + str(st) + ' (Oxford)'
    if species == 'abaumannii_2':
        st = 'PubMLST ST' + str(st) + ' (Pasteur)'
    else:
        st = 'PubMLST ST' + str(st)
    if '-' in st:
        st = 'NA'
    mlst.append(f'{id}\\t{st}\\n')

  with open('mlst_formatted.tsv','w') as outFile:
    outFile.write('Sample\\tMLST Scheme\\n')
    for scheme in mlst:
      outFile.write(scheme)

  """
}

//Kraken Step 1: Run Kraken
process kraken {
  tag "$name"
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: "*_kraken2_report.txt*"

  input:
  set val(name), file(reads) from read_files_kraken

  output:
  tuple name, file("${name}_kraken2_report.txt") into kraken_files, kraken_multiqc

  script:
  """
  kraken2 --db /kraken2-db/minikraken2_v1_8GB --threads ${task.cpus} --report ${name}_kraken2_report.txt --paired ${reads[0]} ${reads[1]}
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

  files = glob.glob("*kraken2_report*")

  results = []
  for file in files:
      sample_id = os.path.basename(file).split(".")[0].replace("_kraken2_report","")
      vals = []
      vals_concat = []
      with open(file,"r") as inFile:
          for line in inFile:
              line = line.strip()
              sline = line.split("\\t")
              if sline[5] == "unclassified":
                  vals.append(sline)
              if sline[3] == "S":
                  vals.append(sline)

      vals_df = DataFrame(vals, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])
      vals_df['Name'] = vals_df['Name'].str.lstrip()
      vals_df = vals_df.sort_values(by=['Percentage'], ascending=False)
      unclass = vals_df[vals_df["Name"]=="unclassified"]
      vals_df = vals_df[vals_df["Name"]!="unclassified"]
      vals_df = vals_df.head(2)
      if len(vals_df) == 0:
          vals_df = vals_df.append(pd.Series(), ignore_index=True)
          vals_df = vals_df.append(pd.Series(), ignore_index=True)
      if len(vals_df) == 1:
          vals_df = vals_df.append(pd.Series(), ignore_index=True)

      vals_concat = pd.concat([unclass,vals_df])
      vals_concat = vals_concat.assign(Sample=sample_id)
      vals_concat = vals_concat[['Sample','Percentage','Name']]
      vals_concat = vals_concat.reset_index(drop=True)

      unclassified = vals_concat.iloc[0]['Percentage'] + "%"
      if str(vals_concat.iloc[1]['Name']) == "nan":
          first_species = "NA"
      else:
          first_species = vals_concat.iloc[1]['Name'] + " (" + vals_concat.iloc[1]['Percentage'] + "%)"
      if str(vals_concat.iloc[2]['Name']) == "nan":
          second_species = "NA"
      else:
          second_species = vals_concat.iloc[2]['Name'] + " (" + vals_concat.iloc[2]['Percentage'] + "%)"

      combined = [[sample_id, unclassified, first_species, second_species]]
      result = DataFrame(combined, columns=['Sample','Unclassified Reads (%)','Primary Species (%)','Secondary Species (%)'])
      results.append(result)

  vals_concat = pd.concat(results)
  vals_concat.to_csv(f'kraken_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//AR Step 1: Run amrfinder
process amrfinder {
  publishDir "${params.outdir}/amrfinder",mode:'copy'

  input:
  file(kraken) from kraken_amr
  file(genomes) from assembled_genomes_ar.collect()

  output:
  file("*_amr.tsv*") into ar_predictions

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import os
  import sys
  import subprocess as sub
  import shlex

  species = ["Acinetobacter_baumannii","Enterococcus_faecalis","Enterococcus_faecium","Staphylococcus_aureus","Staphylococcus_pseudintermedius","Streptococcus_agalactiae","Streptococcus_pneumoniae","Streptococcus_pyogenes"]

  genus = ["Campylobacter","Escherichia","Klebsiella","Salmonella"]

  with open("kraken_results.tsv", "r") as inFile:
      next(inFile)
      for line in inFile:
          line = line.strip().split("\\t")
          sid = line[0]
          taxa = line[2].split(" ")
          taxa_species = taxa[0] + "_" + taxa[1]
          taxa_genus = taxa[0]
          if any(x in taxa_species for x in species):
              outFile = open(f"{sid}_amr.tsv","w")
              cmd = shlex.split(f"amrfinder -n {sid}.contigs.fa --organism {taxa_species}")
              print(cmd)
              sub.Popen(cmd, stdout=outFile).wait()
          if any(x in taxa_genus for x in genus):
              outFile = open(f"{sid}_amr.tsv","w")
              cmd = shlex.split(f"amrfinder -n {sid}.contigs.fa --organism {taxa_genus}")
              print(cmd)
              sub.Popen(cmd, stdout=outFile).wait()
          if taxa_genus == "Shigella":
              outFile = open(f"{sid}_amr.tsv","w")
              cmd = shlex.split(f"amrfinder -n {sid}.contigs.fa --organism Escherichia")
              print(cmd)
              sub.Popen(cmd, stdout=outFile).wait()
          else:
              outFile = open(f"{sid}_amr.tsv","w")
              cmd = shlex.split(f"amrfinder -n {sid}.contigs.fa")
              print(cmd)
              sub.Popen(cmd, stdout=outFile).wait()
  """
}

//AR Step 2: Summarize amrfinder+ results
process amrfinder_summary {
  publishDir "${params.outdir}/amrfinder",mode:'copy'

  input:
  file(predictions) from ar_predictions.collect()

  output:
  file("ar_predictions.tsv")
  file("ar_summary.tsv") into ar_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd

  files = glob.glob("*.tsv")
  dfs = []
  semi_dfs = []
  for file in files:
      sample_id = os.path.basename(file).split(".")[0]
      df = pd.read_csv(file, header=0, delimiter="\\t")
      df.columns=df.columns.str.replace(' ', '_')
      df = df.assign(Sample=sample_id)
      df = df[['Sample','Gene_symbol','%_Coverage_of_reference_sequence','%_Identity_to_reference_sequence']]
      df = df.rename(columns={'%_Coverage_of_reference_sequence':'Coverage','%_Identity_to_reference_sequence':'Identity','Gene_symbol':'Gene'})
      dfs.append(df)
      sample = sample_id
      gene = df['Gene'].tolist()
      gene = ';'.join(gene)
      coverage = df['Coverage'].tolist()
      coverage = ';'.join(map(str, coverage))
      identity = df['Identity'].tolist()
      identity = ';'.join(map(str, identity))
      data = [[sample,gene,coverage,identity]]
      semi_df = pd.DataFrame(data, columns = ['Sample', 'Gene', 'Coverage', 'Identity'])
      semi_dfs.append(semi_df)
  concat_dfs = pd.concat(dfs)
  concat_dfs['Sample'] = concat_dfs['Sample'].str.replace('_amr', '')
  concat_dfs.to_csv('ar_predictions.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  concat_semi_dfs = pd.concat(semi_dfs)
  concat_semi_dfs['Sample'] = concat_semi_dfs['Sample'].str.replace('_amr', '')
  concat_semi_dfs.to_csv('ar_summary.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

process bbduk_summary {
  publishDir "${params.outdir}/bbduk",mode:'copy'

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

  output:
  file('spriggan_report.txt')

  script:
  """
  #!/usr/bin/env python3

  import os
  import glob
  import pandas as pd
  from functools import reduce

  files = glob.glob("*.tsv")

  dfs = []

  for file in files:
      df = pd.read_csv(file, header=0, delimiter='\\t')
      dfs.append(df)

  merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],
                                              how='left'), dfs)

  merged.to_csv('spriggan_report.txt', index=False, sep='\\t', encoding='utf-8')
  """
}

Channel
  .from("$baseDir/multiqc_config.yaml")
  .set { multiqc_config }

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
