# Spriggan
![Spriggan](https://github.com/wslh-bio/spriggan/actions/workflows/spriggan_build.yml/badge.svg)
![GPL-3.0](https://img.shields.io/github/license/wslh-bio/spriggan)
![GitHub Release](https://img.shields.io/github/release/wslh-bio/spriggan)

Spriggan is a [Nextflow](https://www.nextflow.io/) pipeline for the assembly of bacterial whole genome sequence data and identification of antibiotic resistance genes.

### Table of Contents:
[Usage](#usage)  
[Input](#input)  
[Parameters](#parameters)  
[Workflow outline](#workflow-outline)  
[Read trimming and quality assessment](#read-trimming-and-quality-assessment)  
[Genome assembly](#genome-assembly)  
[Assembly quality assessment](#assembly-quality-assessment)  
[Genome coverage](#genome-coverage)  
[Antimicrobial resistance gene detection](#antimicrobial-resistance-gene-detection)  
[MLST scheme](#mlst-scheme)  
[Contamination detection](#contamination-detection)  
[Summary](#summary)  
[Output](#output)  

### Usage
The pipeline is designed to start from raw, paired-end Illumina reads. Start the pipeline using:
```
nextflow spriggan/main.nf --input [path-to-samplesheet] --outdir [path-to-outdir] -profile [docker,singularity,aws]
```

You can specify a version of the pipeline and run it directly from the github repository by using:
```
nextflow wslh-bio/spriggan -r [version] --input [path-to-samplesheet] --outdir [path-to-outdir] -profile [docker,singularity,aws]
```

You can also test the pipeline with example data using `-profile test` or `-profile test_full`:
```
nextflow spriggan/main.nf --outdir [path-to-outdir] -profile test[_full],[docker/singularity]
```
### Input

Spriggan's inputs are paired Illumina FASTQ files for each sample and a comma separated sample sheet containing the sample name, the path to the forward reads file, and the path to the reverse reads file for each sample. An example of the sample sheet's format can be seen in the table below and found [here](https://github.com/wslh-bio/spriggan/blob/main/samplesheets/test_full.csv).

| sample  | fastq_1 | fastq_2 |
| ------------- | ------------- | ------------- |
| sample_name  | /path/to/sample_name_R1.fastq.gz | /path/to/sample_name_R2.fastq.gz |

### Parameters
Spriggan's main parameters and their defaults are shown in the table below:

| Parameter  | Parameter description and default |
| ------------- | ------------- |
| input  | Path to comma-separated file containing information about the samples in the experiment. |
| outdir  | Output directory where the results will be saved. Absolute path must be used for storage on cloud infrastructure. |
| qualitytrimscore | Minimum read quality for trimming (default: 10) |
| trimdirection | Read trimming direction (default: "lr") |
| minlength | Minimum read length for trimming (default: 10) |
| contaminants | Path to FASTA file of contaminating sequences for trimming |
| mincoverage | Minimum coverage threshold to pass a sample (default: 40) |
| kraken_db | Path to Kraken database for classification |
| plus | Use AMRFinderPlus' --plus option (default: false) |
| selected_genes | Genes of interest to pull from AMRFinderPlus output (default: 'NDM\|OXA\|KPC\|IMP\|VIM') |

### Workflow outline

<img src ='/assets/Spriggan.png'>

#### Read trimming and quality assessment
Read trimming and cleaning is performed using [BBtools v38.76](https://jgi.doe.gov/data-and-tools/bbtools/) to trim reads of low quality bases and remove PhiX contamination. Then [FastQC v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used assess the quality of the raw and cleaned reads.

#### Genome assembly
Assembly of the cleaned and trimmed reads is performed using [Shovill v1.1.0](https://github.com/tseemann/shovill).

#### Assembly quality assessment
Quality assessment of the assemblies is performed using [QUAST v5.0.2](http://bioinf.spbau.ru/quast)

#### Genome coverage
Mean and median genome coverage is determined by mapping the cleaned reads back their the assembly using [BWA v0.7.17-r1188](http://bio-bwa.sourceforge.net/) and calculating depth using [samtools v1.10](http://www.htslib.org/)

#### Antimicrobial resistance gene detection
Antimicrobial resistance genes, as well as point mutations, are identified using [AMRFinderPlus v3.10.30](https://github.com/ncbi/amr). Using the plus parameter provides results from the AMRFinderPlus "--plus" option, which includes genes such as virulence factors, stress-response, etc.  
  
Spriggan can also generate a table of results for genes of interest with the selected_genes parameter. Spriggan will search for matches to the gene(s) of interest in the AMRFinderPlus results and make a separate table called **selected_ar_genes.tsv**. The list of genes must be separated by | and enclosed in single quotes in the config file. 

#### MLST scheme
MLST scheme is classified using [MLST v2.17.6](https://github.com/tseemann/mlst). Multiple schemes are available for specific organisms, and STs from all available schemes are reported for those organisms.

#### Contamination detection
Contamination is detected by classifying reads using [Kraken2 v2.0.8](https://ccb.jhu.edu/software/kraken2/) with the Minikraken2_v1_8GB database. A custom Kraken database can be used with the kraken_db parameter.

#### Summary
Results are summarized using [MultiQC v1.11](https://multiqc.info/) and [Pandas v1.3.2](https://pandas.pydata.org/). The main outputs of Spriggan are a csv file named **spriggan_report.csv** and an HTML report file named **spriggan_multiqc_report.html**. The **spriggan_report.csv** file summarizes the results of the QC, classification, and MLST steps of the pipeline. The **spriggan_multiqc_report.html** file contains tables and figures of quality metrics from the FastQC, BBDuk, Samtools, Kraken, and QUAST steps of the pipeline. 

### Output
An example of Spriggan's output directory structure and its output files can be seen below:
```
spriggan_results
├── amrfinder
│   ├── *.amr.tsv
│   ├── *.fa
│   ├── amrfinder_predictions.tsv
│   ├── amrfinder_summary.tsv
│   └── selected_ar_genes.tsv
├── bbduk
│   ├── *.fastq.gz
│   ├── *.adapter.stats.txt
│   ├── *.bbduk.log
│   ├── *.trim.txt
│   └── bbduk_results.tsv
├── coverage
│   └── coverage_stats.tsv
├── fastqc
│   ├── *.html
│   ├── *.zip
│   └── fastqc_summary.tsv
├── kraken
│   ├── *.kraken2.txt
│   ├── kraken_results.tsv
│   └── kraken2.log
├── mlst
│   ├── *.alleles.tsv
│   ├── *.mlst.tsv
│   └── mlst_results.tsv
├── multiqc
│   ├── multiqc_data
│   │   ├── *.json
│   │   ├── *.txt
│   │   └── multiqc.log
│   ├── multiqc_plots
│   │   ├── pdf
│   │   │   └── *.pdf
│   │   ├── png
│   │   │   └── *.png
│   │   └── svg
│   │       └── *.svg
│   └── spriggan_multiqc_report.html
├── pipeline_info
│   ├── *.html
│   ├── *.txt
│   ├── samplesheet.valid.csv
│   └── software_versions.yml
├── quast
│   ├── *.quast.report.tsv
│   ├── *.transposed.quast.report.tsv
│   └── quast_results.tsv
├── results
│   └── spriggan_report.csv
├── samtools
│   ├── *.bam
│   ├── *.depth.tsv
│   └── *.stats.txt
└── shovill
    ├── *.contigs.fa
    ├── *.sam
    └── shovill_output
          ├── contigs.gfa
          ├── shovill.corrections
          ├── shovill.log
          └── spades.fasta
```
Notable output files:  
**spriggan_report.csv** - Summary table of each step in Spriggan  
**spriggan_multiqc_report.html** - HTML report generated by MultiQC  
**\*.contigs.fa** - Shovill assembly for each sample  
**\*.amr.tsv** - AMR genes identified in each sample by AMRFinderPlus  
**\*.mlst.tsv** - MLST scheme identified for each sample

### Authors
[Kelsey Florek](https://github.com/k-florek), WSLH Senior Genomics and Data Scientist  
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatician and Data Scientist

### Citations
```
K. Florek & A.C. Shockey (2021). Spriggan (Version 1.5.3) [https://github.com/wslh-bio/spriggan].
```
