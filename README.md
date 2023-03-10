# Spriggan
![Spriggan](https://github.com/wslh-bio/spriggan/actions/workflows/spriggan_build.yml/badge.svg)
![GPL-3.0](https://img.shields.io/github/license/wslh-bio/spriggan)
![GitHub Release](https://img.shields.io/github/release/wslh-bio/spriggan)

Spriggan is a [NextFlow](https://www.nextflow.io/) pipeline used for assembly of bacterial whole genome sequence data and identification of antibiotic resistance genes.

### Table of Contents:
[Usage](#usage)  
[Workflow outline](#workflow-outline)  
[Read trimming and quality assessment](#read-trimming-and-quality-assessment)  
[Genome assembly](#genome-assembly)  
[Assembly quality assessment](#assembly-quality-assessment)  
[Genome coverage](#genome-coverage)  
[Antimicrobial resistance gene detection](#antimicrobial-resistance-gene-detection)  
[MLST scheme](#mlst-scheme)  
[Contamination detection](#contamination-detection)  
[Summary](#summary)  
[Output files](#output-files)  

### Usage
The pipeline is designed to start from raw Illumina reads. All reads must be in the same directory. Then start the pipeline using:
```
nextflow main.nf --reads [path-to-reads]
```

You can also test the pipeline with example data using `--test`, note this requires NextFlow version `21.07.0-edge` or greater:
```
nextflow main.nf --test
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
Antimicrobial resistance genes, as well as point mutations, are identified using [AMRFinderPlus v3.10.30](https://github.com/ncbi/amr). Spriggan can generate a table of AMRFinderPlus results for AR genes of interest with the selected_genes parameter. Spriggan will search for matches to the AR genes of interest in the AMRFinderPlus results and make a table called 'selected_ar_genes.tsv.' The list of genes must be separated by | and enclosed in single quotes in the config file. By default the selected_genes parameter is set to: 'NDM|OXA|KPC|IMP|VIM'

#### MLST scheme
MLST scheme is classified using [MLST v2.17.6](https://github.com/tseemann/mlst). Multiple schemes are available for specific organisms, and STs from all available schemes are reported for those organisms.

#### Contamination detection
Contamination is detected by classifying reads using [Kraken2 v2.0.8](https://ccb.jhu.edu/software/kraken2/) with the Minikraken2_v1_8GB database.

#### Summary
Results are summarized using [MultiQC v1.11](https://multiqc.info/) and [Pandas v1.3.2](https://pandas.pydata.org/).

### Output files

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

**amrfinder_predictions.tsv** - Long-format table of AMRFinderPlus results  
**amrfinder_summary.tsv** - Wide-format table of AMRFinderPlus results  
**\*.amr.tsv** - Raw AMRFinderPlus output for each sample  
**selected_ar_genes.tsv** - Table of AMRFinderPlus results for genes selected by user  
**\*.contigs.fa** - Shovill assembly for each sample  
**fastqc_summary.tsv** - Summary table of FastQC results  
**\*.html** - HTML files of FastQC results  
**\*.zip** - Zipped folders of FastQC output  
**kraken_results.tsv** - Summary table of Kraken results  
**\*.kraken2.txt** - Report of Kraken results for each sample  
**\*.bam** - Alignments to an assembly BAM format  
**\*.bai** - Index file of alignments to an assembly  
**coverage_stats.tsv** - Summary table of mean and median coverage calculated with Samtools depth  
**\*.depth.tsv** - Raw Samtools depth output for reads mapped to an assembly  
**\*.sam** - Alignments to an assembly SAM format  
**\*.stats.txt** - Output of Samtools stats for reads mapped to an assembly  
**\*.alleles.tsv** - Raw MLST alleles for each sample  
**mlst_results.tsv** - Summary table of MLST results  
**\*.mlst.tsv** - Raw MLST output for each sample  
**quast_results.tsv** - Summary table of QUAST results  
**\*.quast.report.tsv** - Long-format QUAST results for each sample  
**\*.transposed.quast.report.tsv** - Wide-format QUAST results for each sample  
**multiqc_report.html** - HTML report generated by MultiQC  
**spriggan_report.csv** - Summary table of each step in Spriggan  
**bbduk_results.tsv** - Summary table of trimming with BBDuk  
**\*\_clean\_\*** - Trimmed and cleaned reads  
**\*.trim.txt** - Trimming results from BBDuk each sample  

### Authors
[Kelsey Florek](https://github.com/k-florek), WSLH Senior Genomics and Data Scientist  
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatician and Data Scientist
