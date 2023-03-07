/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BBDUK                         } from '../modules/local/bbduk.nf'
include { BBDUK_SUMMARY                 } from '../modules/local/bbduk_summary.nf'
include { FASTQC                        } from '../modules/local/fastqc.nf'
include { FASTQC_SUMMARY                } from '../modules/local/fastqc_summary.nf'
include { SHOVILL                       } from '../modules/local/shovill.nf'
include { SAMTOOLS                      } from '../modules/local/samtools.nf'
include { COVERAGE_STATS                } from '../modules/local/coverage_stats.nf'
include { QUAST                         } from '../modules/local/quast.nf'
include { QUAST_SUMMARY                 } from '../modules/local/quast_summary.nf'
include { MLST                          } from '../modules/local/mlst.nf'
include { MLST_SUMMARY                  } from '../modules/local/mlst_summary.nf'
include { KRAKEN                        } from '../modules/local/kraken.nf'
include { KRAKEN_SUMMARY                } from '../modules/local/kraken_summary.nf'
include { AMRFINDER_SETUP               } from '../modules/local/amrfinder_setup.nf'
include { AMRFINDER                     } from '../modules/local/amrfinder.nf'
include { AMRFINDER_SUMMARY             } from '../modules/local/amrfinder_summary.nf'
include { RESULTS                       } from '../modules/local/results.nf'
include { MULTIQC                       } from '../modules/local/multiqc.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SPRIGGAN {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    INPUT_CHECK.out.reads
        .branch {
            ntc: it[0]['id'].contains('NTC')
            sample: !it[0]['id'].contains('NTC')
        }
        .set{ ch_input_reads }

    //
    // MODULE: BBDUK
    //
    BBDUK (
        ch_input_reads.sample,
        params.contaminants
    )
    ch_versions = ch_versions.mix(BBDUK.out.versions.first())

    //
    // MODULE: BBDUK_SUMMARY
    //
    BBDUK_SUMMARY (
        BBDUK.out.bbduk_trim.collect()
    )

    //
    // MODULE: FASTQC
    //
    FASTQC (
        BBDUK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: FASTQC_SUMMARY
    //
    FASTQC_SUMMARY (
        FASTQC.out.zip.collect{it[1]}
    )

    //
    // MODULE: SHOVILL
    //
    SHOVILL (
        BBDUK.out.reads
    )
    ch_versions = ch_versions.mix(SHOVILL.out.versions.first())

    //
    // MODULE: SAMTOOLS
    //
    SAMTOOLS (
        SHOVILL.out.sam_files
    )
    ch_versions = ch_versions.mix(SAMTOOLS.out.versions.first())

    //
    // MODULE: COVERAGE_STATS
    //
    COVERAGE_STATS (
        SAMTOOLS.out.cov_files.collect()
    )

    //
    // MODULE: QUAST
    //
    QUAST (
        SHOVILL.out.contigs
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    //
    // MODULE: QUAST_SUMMARY
    //
    QUAST_SUMMARY (
        QUAST.out.transposed_report.collect()
    )

    //
    // MODULE: MLST
    //
    MLST (
        SHOVILL.out.contigs
    )
    ch_versions = ch_versions.mix(MLST.out.versions.first())

    //
    // MODULE: MLST_SUMMARY
    //
    MLST_SUMMARY (
        MLST.out.mlst_files.collect()
    )

    //
    // MODULE: KRAKEN
    //

    if (params.kraken_db != null) {
        Channel
            .fromPath(params.kraken_db)
            .set { kraken_db }
    } else {
        kraken_db = file("$baseDir/assets/empty.txt",checkIfExists:true)
    }

    KRAKEN (
        BBDUK.out.reads,
        kraken_db
    )
    ch_versions = ch_versions.mix(KRAKEN.out.versions.first())

    //
    // MODULE: KRAKEN_SUMMARY
    //
    KRAKEN_SUMMARY (
        KRAKEN.out.kraken_results.collect()
    )

    //
    // MODULE: AMRFINDER_SETUP
    //
    AMRFINDER_SETUP (
        KRAKEN_SUMMARY.out.kraken_tsv,
        SHOVILL.out.contigs
    )

    //
    // MODULE: AMRFINDER
    //
    AMRFINDER (
        AMRFINDER_SETUP.out.amrfinder_input
    )
    ch_versions = ch_versions.mix(AMRFINDER.out.versions.first())

    //
    // MODULE: AMRFINDER_SUMMARY
    //
    AMRFINDER_SUMMARY (
        AMRFINDER.out.amrfinder_predictions.collect()
    )

    //
    // MODULE: RESULTS
    //
    RESULTS (
        BBDUK_SUMMARY.out.bbduk_tsv,
        COVERAGE_STATS.out.coverage_tsv,
        QUAST_SUMMARY.out.quast_tsv,
        MLST_SUMMARY.out.mlst_tsv,
        KRAKEN_SUMMARY.out.kraken_tsv,
        AMRFINDER_SUMMARY.out.amrfinder_tsv,
        AMRFINDER_SUMMARY.out.selected_ar_tsv,
        KRAKEN.out.versions.first(),
        AMRFINDER.out.versions.first()
    )

    //
    // MODULE: WORKFLOW_TEST
    //
    /*
    ch_valid_dataset = Channel.fromPath("$projectDir/test-dataset/validation/spntypeid_report_valid.csv", checkIfExists: true)
    WORKFLOW_TEST (
        ch_valid_dataset.collect(),
        RESULTS.out.result_csv
    )
    */

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSpriggan.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSpriggan.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BBDUK.out.bbduk_adapters.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BBDUK.out.bbduk_trim.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS.out.stats_multiqc.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN.out.kraken_results.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.result.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
