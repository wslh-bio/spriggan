# Spriggan

## Note: Spriggan is currently a work in progress and may not be fully functional

Spriggan is a [NextFlow](https://www.nextflow.io/) pipeline that assembles bacterial whole genome sequence data and identifies antibiotic resistance genes.

### Table of Contents:
[Usage](#using-the-pipeline)  
[Workflow outline](#workflow-outline)  
[Read trimming and quality assessment](#read-trimming-and-quality-assessment)  
[Genome assembly](#genome-assembly)  
[Assembly quality assessment](#assembly-quality-assessment)  
[Genome coverage](#genome-coverage)  
[MLST scheme](#mlst-scheme)  
[Contamination detection](#contamination-detection)
[Output](#output-files)  
[Dependencies](#dependencies)   

### Using the pipeline
The pipeline is designed to start from raw Illumina short reads. All reads must be in the same directory. Then start the pipeline using `nextflow sprriggan.nf --reads [path-to-reads]`.

### Workflow outline

#### Read trimming and quality assessment
Read trimming and cleaning is performed using [BBtools v38.76](https://jgi.doe.gov/data-and-tools/bbtools/) to trim reads of low quality bases and remove PhiX contamination. Then [FastQC v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is used assess the quality of the raw and cleaned reads.

#### Genome assembly
Assembly of the cleaned and trimmed reads is performed using [Shovill v1.1.0](https://github.com/tseemann/shovill).

#### Assembly quality assessment
Quality assessment of the assemblies is performed using [QUAST v5.0.2](http://bioinf.spbau.ru/quast)

#### Genome coverage
Mean and median genome coverage is determined by mapping the cleaned reads back their the assembly using [BWA v0.7.17-r1188](http://bio-bwa.sourceforge.net/) and calculating depth using [samtools v1.10](http://www.htslib.org/)

#### Antimicrobial resistance gene detection
Antimicrobial resistance genes, as well as point mutations, are identified using [AMRFinderPlus v3.1.1](https://github.com/ncbi/amr).

#### MLST scheme
MLST scheme is classified using [MLST v2.17.6](https://github.com/tseemann/mlst). Multiple schemes are available for specific organisms, and STs from all available schemes are reported for those organisms.

#### Contamination detection
Contamination is detected by classifying reads using [Kraken2 v2.0.8](https://ccb.jhu.edu/software/kraken2/) with the Minikraken database.

### Output files

```
spriggan_results
├── alignments
│   ├── *.sam
│   └── *.bam
├── amrfinder
│   ├── ar_predictions.tsv
│   ├── ar_summary.tsv
│   └── *amr.tsv
├── assembled
│   └── *.contigs.fa
├── bbduk
│   └── bbduk_results.tsv
├── coverage
│   ├── coverage_stats.tsv
│   └── *depth.tsv
├── fastqc
│   ├── fq_summary.txt
│   └── *.html
├── kraken
│   ├── kraken_results.tsv
│   └── *kraken2_report.txt
├── mlst
│   ├── mlst_results.tsv
│   └── *.tsv
├── quast
│   ├── quast_results.tsv
│   └── *.quast.tsv
├── multiqc_report.html
├── spriggan_report.txt
└── trimming
    └── *.trim.txt
```

### Authors
[Kelsey Florek](https://github.com/k-florek), WSLH Bioinformatics Scientist  
[Abigail Shockey](https://github.com/AbigailShockey), WSLH Bioinformatics Scientist
