# Spriggan
![Spriggan](https://github.com/wslh-bio/spriggan/actions/workflows/spriggan_build.yml/badge.svg)
![GPL-3.0](https://img.shields.io/github/license/wslh-bio/spriggan)
![GitHub Release](https://img.shields.io/github/release/wslh-bio/spriggan)

Spriggan is a [NextFlow](https://www.nextflow.io/) pipeline used for assembly of bacterial whole genome sequence data and identification of antibiotic resistance genes.

### Table of Contents:
[Usage](#using-the-pipeline)  
[Workflow outline](#workflow-outline)  
[Read trimming and quality assessment](#read-trimming-and-quality-assessment)  
[Genome assembly](#genome-assembly)  
[Assembly quality assessment](#assembly-quality-assessment)  
[Genome coverage](#genome-coverage)  
[Antimicrobial resistance gene detection](#antimicrobial-resistance-gene-detection)  
[MLST scheme](#mlst-scheme)  
[Contamination detection](#contamination-detection)                                                                                                                                   
[Output](#output-files)  

### Using the workflow
The pipeline is designed to start from raw Illumina reads. All reads must be in the same directory. Then start the pipeline using:
```
nextflow sprriggan.nf --reads [path-to-reads]
```

You can also test the pipeline with example data using `--test`, note this requires NextFlow version `21.07.0-edge` or greater:
```
nextflow sprriggan.nf --test
```

You can specify a version of the pipeline and run it directly from the github repository by using:
```
nextflow wslh-bio/spriggan -r <version> --reads [path-to-reads]
```

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
Antimicrobial resistance genes, as well as point mutations, are identified using [AMRFinderPlus v3.1.1](https://github.com/ncbi/amr). Spriggan can generate a table of AMRFinderPlus results for AR genes of interest with the selected_genes parameter. Spriggan will search for matches to the AR genes of interest in the AMRFinderPlus results and make a table called 'selected_ar_genes.tsv.' The list of genes must be separated by | and enclosed in single quotes in the config file. By default the selected_genes parameter is set to: 'NDM|OXA|KPC|IMP|VIM'

#### MLST scheme
MLST scheme is classified using [MLST v2.17.6](https://github.com/tseemann/mlst). Multiple schemes are available for specific organisms, and STs from all available schemes are reported for those organisms.

#### Contamination detection
Contamination is detected by classifying reads using [Kraken2 v2.0.8](https://ccb.jhu.edu/software/kraken2/) with the Minikraken2_v1_8GB database.

### Output files

```
spriggan_results
????????? alignments
???   ????????? *.(sorted).bam(.bai)
???   ????????? *.sam
???   ????????? *.txt
????????? amrfinder
???   ????????? ar_predictions.tsv
???   ????????? ar_summary.tsv
???   ????????? *.amr.tsv
???   ????????? selected_ar_genes.tsv
????????? assembled
???   ????????? *.contigs.fa
????????? coverage
???   ????????? coverage_stats.tsv
???   ????????? *.depth.tsv
????????? fastqc
???   ????????? fastqc_summary.tsv
???   ????????? *.html
???   ????????? zips
???       ????????? *.zip
????????? kraken
???   ????????? kraken_results.tsv
???   ????????? *.kraken2.txt
????????? mlst
???   ????????? mlst_results.tsv
???   ????????? *.mlst.tsv
????????? quast
???   ????????? quast_results.tsv
???   ????????? *.quast.tsv
????????? multiqc_report.html
????????? spriggan_report.csv
????????? trimming
    ????????? bbduk_results.tsv
    ????????? *.trim.txt
```

**\*.(sorted).bam(.bai)** - Alignments in BAM format (sorted and unsorted; includes \*.bai index files)  
**\*.sam** - Alignments in SAM format  
**\*.txt** - Samtools stats output  
**ar_predictions.tsv** - Long-format table of AMRFinderPlus results  
**ar_summary.tsv** - Wide-format table of AMRFinderPlus results  
**\*.amr.tsv** - Raw AMRFinderPlus output for each sample  
**selected_ar_genes.tsv** - Table of AMRFinderPlus results for genes selected by user  
**\*.contigs.fa** - Shovill assembly for each sample  
**coverage_stats.tsv** - Summary table of mean and median coverage calculated with Samtools depth  
**\*.depth.tsv** - Raw Samtools depth output for each sample  
**fastqc_summary.tsv** - Summary table of FastQC results  
**\*.html** - HTML files of FastQC results  
**\*.zip** - Zipped folders of FastQC output    
**kraken_results.tsv** - Summary table of Kraken results  
**\*.kraken2.txt** - Report of Kraken results for each sample  
**mlst_results.tsv** - Summary table of MLST results  
**\*.mlst.tsv** - Raw MLST output for each sample  
**quast_results.tsv** - Summary table of QUAST results  
**\*.quast.tsv** - QUAST results for each sample  
**multiqc_report.html** - HTML report generated by MultiQC  
**spriggan_report.csv** - Summary table of each step in Spriggan  
**bbduk_results.tsv** - Summary table of trimming with BBDuk  
**\*.trim.txt** - Trimming results from BBduk each sample

### Authors
[Kelsey Florek](https://github.com/k-florek), WSLH Senior Genomics and Data Scientist  
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatics and Data Scientist
