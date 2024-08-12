/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//Channel.value(params).view()
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.filter_table,
    params.binddir
]
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
========================================================================================
    INPUT AND VARIABLES
========================================================================================
*/

ch_versions = Channel.empty()
ch_fasta = Channel.empty()
ch_otu_table = Channel.empty()
ch_lca_input_table = Channel.empty()
ch_curated_table = Channel.empty()
ch_curated_fasta = Channel.empty()
ch_pngs = Channel.empty()
ch_raw_data = Channel.empty()
ch_ulimit = Channel.empty()
ch_db = Channel.fromPath(params.dbfiles).collect()
ch_mitodb = Channel.fromPath(params.mitodbfiles).collect()
ch_caabmap = Channel.fromPath(params.caabmap).collect()

if (params.filter_table) {
    ch_filter = file(params.filter_table, checkIfExists: true)
} else {
    ch_filter = []
}

if (!params.skip_demux) {
    if (params.raw_data) {
        ch_raw_data = Channel.fromFilePairs(params.raw_data, checkIfExists: true)
    } else {
        exit 1, 'Raw data not specified!'
    }

    if (params.ulimit) {
        ch_ulimit = Channel.of(params.ulimit)
    } else {
        exit 1, 'Ulimit not specified!'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BLAST_BLASTN                } from '../modules/local/blast/blastn/main.nf'
include { CONCAT_BLASTN_RESULTS       } from '../modules/local/custom/concat_blastn_results/main.nf'
include { LCA                         } from '../modules/local/lca/main.nf'
include { PHYLOSEQ                    } from '../modules/local/phyloseq/main.nf'
include { REMOVE_DUPS                 } from '../modules/local/custom/removedups/main.nf'
include { OCOMNBC                     } from '../modules/local/custom/ocomnbc/main.nf'
include { MARKDOWN_REPORT             } from '../modules/local/custom/markdownreport/main.nf'
include { SEQTK_TRIM                  } from '../modules/local/seqtk/trim/main.nf'
include { CUTADAPT as CUTADAPT_TRIM   } from '../modules/local/cutadapt/main.nf'
include { SEQKIT_STATS as FINAL_STATS } from '../modules/local/seqkit_stats/main.nf'
include { OBITOOLS3_WORKFLOW          } from '../subworkflows/local/obitools3_workflow'
include { CUTADAPT_WORKFLOW           } from '../subworkflows/local/cutadapt_workflow'
include { ASV_WORKFLOW                } from '../subworkflows/local/asv_workflow'
include { INPUT_CHECK                 } from '../subworkflows/local/input_check'
include { LULU_WORKFLOW               } from '../subworkflows/local/lulu_workflow'
include { ZOTU_WORKFLOW               } from '../subworkflows/local/zotu_workflow'
include { POSTDEMUX_WORKFLOW          } from '../subworkflows/local/postdemux_workflow'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { DOWNLOAD_AQUAMAPS           } from '../modules/local/custom/download_aquamaps/main'
include { GET_AQUAMAP_PROBS           } from '../modules/local/custom/getaquamapprobs/main'
include { GET_CAAB_PROBS              } from '../modules/local/custom/getcaabprobs/main'
include { PRIMER_CONTAM_STATS         } from '../modules/local/custom/primercontamstats/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow OCEANOMICS_AMPLICON {

    if (params.assay_readlength) {
        params.seqtk_trim       = params.assay_params[params.assay_readlength]["seqtk_trim"]
        params.seqtk_length     = params.assay_params[params.assay_readlength]["seqtk_length"]
        params.asv_min_overlap  = params.assay_params[params.assay_readlength]["asv_min_overlap"]
        params.lca_qcov         = params.assay_params[params.assay_readlength]["lca_qcov"]
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Demultiplex with cutadapt or obitools3
    //
    if (!params.skip_demux) {
        if (!params.obi3_demux) {
            CUTADAPT_WORKFLOW (
                ch_input,
                ch_raw_data,
                ch_ulimit
            )
            ch_versions = ch_versions.mix(CUTADAPT_WORKFLOW.out.versions)

            ch_demux_reads = CUTADAPT_WORKFLOW.out.reads
            ch_raw_stats = CUTADAPT_WORKFLOW.out.raw_stats
            ch_raw_stats_collected = ch_raw_stats.map{ it = it[1] }.collect()
            ch_assigned_stats = CUTADAPT_WORKFLOW.out.assigned_stats
            ch_assigned_stats_collected = ch_assigned_stats.map{ it = it[1] }.collect()

        } else {
            OBITOOLS3_WORKFLOW (
                ch_input,
                ch_raw_data
            )
            ch_versions = ch_versions.mix(OBITOOLS3_WORKFLOW.out.versions)

            ch_demux_reads = OBITOOLS3_WORKFLOW.out.reads
            ch_raw_stats = OBITOOLS3_WORKFLOW.out.raw_stats
            ch_raw_stats_collected = ch_raw_stats.map{ it = it[1] }.collect()
            ch_assigned_stats = OBITOOLS3_WORKFLOW.out.raw_stats
            ch_assigned_stats_collected = ch_assigned_stats.map{ it = it[1] }.collect()
        }

        POSTDEMUX_WORKFLOW (
            ch_demux_reads,
            ch_input,
            ch_raw_data
        )
        ch_versions = ch_versions.mix(POSTDEMUX_WORKFLOW.out.versions)

        ch_reads = POSTDEMUX_WORKFLOW.out.reads
        ch_missing = POSTDEMUX_WORKFLOW.out.missing_samples

    } else {
        ch_reads = INPUT_CHECK.out.reads
        ch_raw_stats_collected = []
        ch_assigned_stats_collected = []
        ch_missing = []
    }

    // MODULE: Trim primer sequences
    if (! params.skip_primertrim) {
        CUTADAPT_TRIM (
            ch_reads,
            [],
            [],
            params.ulimit
        )
        ch_versions = ch_versions.mix(CUTADAPT_TRIM.out.versions.first())

        ch_primertrimmed_reads = CUTADAPT_TRIM.out.reads
    } else {
        ch_primertrimmed_reads = ch_reads
    }

    //
    // MODULE: Trim right ends of reads that are too long
    //
    if (params.seqtk_trim) {
        SEQTK_TRIM (
            ch_primertrimmed_reads
        )
        ch_versions = ch_versions.mix(SEQTK_TRIM.out.versions.first())

        ch_trimmed_reads = SEQTK_TRIM.out.reads
    } else {
        ch_trimmed_reads = ch_primertrimmed_reads
    }

    ch_trimmed_reads_collected = ch_trimmed_reads.map{it = it[1]}.collect().map{it = ["concat", it]}

    FINAL_STATS (
        ch_trimmed_reads_collected,
        "final"
    )

    ch_final_stats = FINAL_STATS.out.stats
    ch_final_stats_collected = ch_final_stats.map{ it = it[1] }.collect()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_trimmed_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // SUBWORKFLOW: ASV creation with DADA2
    //
    if (!params.skip_asvs) {
        ASV_WORKFLOW (
            ch_trimmed_reads
        )
        ch_versions = ch_versions.mix(ASV_WORKFLOW.out.versions)

        ch_fasta           = ch_fasta.mix(ASV_WORKFLOW.out.fasta)
        ch_otu_table       = ch_otu_table.mix(ASV_WORKFLOW.out.table)
        ch_lca_input_table = ch_lca_input_table.mix(ASV_WORKFLOW.out.lca_input_table)

        ch_pngs            = ch_pngs.mix(ASV_WORKFLOW.out.quality_raw.randomSample(3, 4).map { it = it[1] })
        ch_pngs            = ch_pngs.mix(ASV_WORKFLOW.out.quality_filt.randomSample(3, 4).map { it = it[1] })
        ch_pngs            = ch_pngs.mix(ASV_WORKFLOW.out.errors_plot)
        ch_pngs            = ch_pngs.mix(ASV_WORKFLOW.out.seq_distribution)
        ch_pngs            = ch_pngs.mix(ASV_WORKFLOW.out.track_reads)
    }

    //
    // SUBWORKFLOW: ZOTU creation with Vsearch, Usearch32, or Usearch64
    //
    if (!params.skip_zotus) {
        ZOTU_WORKFLOW (
            ch_trimmed_reads
        )
        ch_versions = ch_versions.mix(ZOTU_WORKFLOW.out.versions)

        ch_fasta           = ch_fasta.mix(ZOTU_WORKFLOW.out.fasta)
        ch_otu_table       = ch_otu_table.mix(ZOTU_WORKFLOW.out.table)
        ch_lca_input_table = ch_lca_input_table.mix(ZOTU_WORKFLOW.out.lca_input_table)
    }

    //
    // SUBWORKFLOW: LULU workflow
    //
    if (!params.skip_lulu) {
        LULU_WORKFLOW (
            ch_fasta,
            ch_otu_table
        )
        ch_versions = ch_versions.mix(LULU_WORKFLOW.out.versions.first())

        ch_curated_table = ch_curated_table.mix(LULU_WORKFLOW.out.curated_table)
        ch_curated_fasta = ch_curated_fasta.mix(LULU_WORKFLOW.out.curated_fasta)
    } else {
        ch_curated_table = ch_lca_input_table
        ch_curated_fasta = ch_fasta
    }

    ch_curated_fasta_split = ch_curated_fasta
        .splitFasta( by: 1000, file: true )

    //
    // MODULE: Run Blastn
    //
    if (!params.skip_classification) {
        BLAST_BLASTN (
            ch_curated_fasta_split,
            ch_db
        )
        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

        CONCAT_BLASTN_RESULTS (
            BLAST_BLASTN.out.txt.groupTuple()
        )

        ch_lca_input = ch_curated_table.join(CONCAT_BLASTN_RESULTS.out.txt)

        //
        // MODULE: run LCA
        //
        LCA (
            ch_lca_input
        )
        ch_versions = ch_versions.mix(LCA.out.versions.first())
        REMOVE_DUPS (
            LCA.out.lca_output
        )
        ch_versions = ch_versions.mix(REMOVE_DUPS.out.versions.first())

        //
        // MODULE: Naive Bayes Classifier
        //
        OCOMNBC (
            ch_curated_fasta
        )
        ch_versions = ch_versions.mix(OCOMNBC.out.versions.first())

        ch_phyloseq_input = ch_otu_table.join(REMOVE_DUPS.out.tsv.join(OCOMNBC.out.nbc_output))

        //
        // MODULE: Create Phyloseq object
        //
        PHYLOSEQ (
            ch_phyloseq_input,
            ch_input,
            ch_filter
        )
        ch_versions = ch_versions.mix(PHYLOSEQ.out.versions.first())

        //
        // MODULE: Download aquamap nc files
        //
        DOWNLOAD_AQUAMAPS (
            PHYLOSEQ.out.phyloseq_object
        )

        //
        // MODULE: Get the aquamap probabilities for each species in each ASV
        //
        GET_AQUAMAP_PROBS (
            PHYLOSEQ.out.phyloseq_object,
            DOWNLOAD_AQUAMAPS.out.nc_files
        )

        //
        // MODULE: Get the CAAB probabilities for each species in each ASV
        //
        //GET_CAAB_PROBS (
        //    PHYLOSEQ.out.phyloseq_object,
        //    ch_caabmap
        //)

        ch_taxa_collected = PHYLOSEQ.out.final_taxa.collect()
    } else {
        ch_taxa_collected = Channel.empty()
    }

    //
    // MODULE: Collect stats about primer sequences found in ASVs/ZOTUs
    //
    PRIMER_CONTAM_STATS (
        ch_curated_fasta,
        params.fw_primer,
        params.rv_primer
    )

    //
    // MODULE: Create Markdown reports
    //
    MARKDOWN_REPORT (
        ch_final_stats_collected,
        ch_raw_stats_collected,
        ch_assigned_stats_collected,
        ch_taxa_collected.ifEmpty([]),
        ch_pngs.collect(),
        ch_missing,
        ch_input,
        PRIMER_CONTAM_STATS.out.txt.map{ it = it[1] }.collect()
    )
    ch_versions = ch_versions.mix(MARKDOWN_REPORT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowOceanOmicsAmplicon.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowOceanOmicsAmplicon.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MARKDOWN_REPORT.out.html.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
