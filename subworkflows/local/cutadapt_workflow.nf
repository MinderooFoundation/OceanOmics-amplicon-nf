/*
 * Demultiplex with Cutadapt
 */

include { CUTADAPT as CUTADAPT_DEMUX_UC  } from '../../modules/local/cutadapt/main.nf'
include { CUTADAPT as CUTADAPT_DEMUX_UDI } from '../../modules/local/cutadapt/main.nf'
include { CREATE_DEMUX_DEPENDENCIES      } from '../../modules/local/custom/createdemuxdependencies/main.nf'
include { RENAME                         } from '../../modules/local/custom/rename/main.nf'
include { SEQKIT_STATS as \
            ASSIGNED_STATS;
            SEQKIT_STATS as \
            UNKNOWN_STATS;
            SEQKIT_STATS as \
            RAW_STATS                    } from '../../modules/local/seqkit_stats/main.nf'
include { CONCAT                         } from '../../modules/local/custom/concat/main.nf'

workflow CUTADAPT_WORKFLOW {
    take:
    ch_input    // channel: [ samplesheet.csv ]
    ch_raw_data // channel: [ prefix, [ reads ] ]
    ch_ulimit   // channel: [ ulimit ]

    main:
    ch_versions = Channel.empty()

    // MODULE: Check stats of raw data
    RAW_STATS (
        ch_raw_data,
        "raw"
    )
    ch_versions = ch_versions.mix(RAW_STATS.out.versions)

    //
    // MODULE: Create files needed for demultiplexing and renaming of files
    //
    CREATE_DEMUX_DEPENDENCIES (
        ch_input,
        params.demux_udi
    )
    ch_versions = ch_versions.mix(CREATE_DEMUX_DEPENDENCIES.out.versions)

    // MODULE: Demultiplex
    if (! params.demux_udi) {
        CUTADAPT_DEMUX_UC (
            ch_raw_data,
            CREATE_DEMUX_DEPENDENCIES.out.fw_index,
            CREATE_DEMUX_DEPENDENCIES.out.rv_index,
            params.ulimit
        )
        ch_versions = ch_versions.mix(CUTADAPT_DEMUX_UC.out.versions)

        // MODULE: Rename the Cutadapt output files
        RENAME (
            CUTADAPT_DEMUX_UC.out.reads,
            CREATE_DEMUX_DEPENDENCIES.out.sample_rename
        )
        ch_versions = ch_versions.mix(RENAME.out.versions)

        // MODULE: Check stats of reads that couldn't be assigned to samples after demultiplexing
        UNKNOWN_STATS (
            RENAME.out.unknown,
            "unknown"
        )
        ch_versions = ch_versions.mix(UNKNOWN_STATS.out.versions)

        // MODULE: Trim leftover primers and concatenate files so that there is only one R1 and one R2 file per sample
        // NOTE: This process initially was being used to trim primers and concatenate files, now it's just concatenating.
        //       the trimming of primers is now being done by Cutadapt later in the pipeline
        CONCAT (
            RENAME.out.reads,
            ch_input
        )
        ch_versions = ch_versions.mix(CONCAT.out.versions)

        ch_demuxed_reads = CONCAT.out.reads
    } else {
        CUTADAPT_DEMUX_UDI (
            ch_raw_data,
            CREATE_DEMUX_DEPENDENCIES.out.fw_index,
            CREATE_DEMUX_DEPENDENCIES.out.rv_index,
            params.ulimit
        )
        ch_versions = ch_versions.mix(CUTADAPT_DEMUX_UDI.out.versions)

        ch_demuxed_reads = CUTADAPT_DEMUX_UDI.out.reads
    }

     // MODULE: Check stats of reads assigned to samples after demultiplexing
    ASSIGNED_STATS (
        ch_demuxed_reads,
        "assigned"
    )
    ch_versions = ch_versions.mix(ASSIGNED_STATS.out.versions)

    emit:
    reads           = ch_demuxed_reads                    // channel: [ val(meta), reads ]
    assigned_stats  = ASSIGNED_STATS.out.stats            // channel: [ val(meta), stats ]
    //unknown_stats   = UNKNOWN_STATS.out.stats             // channel: [ val(meta), stats ]
    raw_stats       = RAW_STATS.out.stats                 // channel: [ val(meta), stats ]
    versions        = ch_versions                         // channel: [ versions.yml ]
}
