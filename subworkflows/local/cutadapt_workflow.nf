/*
 * Demultiplex with Cutadapt
 */

include { CUTADAPT                  } from '../../modules/local/cutadapt/main.nf'
include { CREATE_DEMUX_DEPENDENCIES } from '../../modules/local/custom/createdemuxdependencies/main.nf'
include { RENAME                    } from '../../modules/local/custom/rename/main.nf'
include { SEQKIT_STATS as \
          ASSIGNED_STATS;
          SEQKIT_STATS as \
          UNKNOWN_STATS;
          SEQKIT_STATS as \
          RAW_STATS;
          SEQKIT_STATS as \
          FINAL_STATS               } from '../../modules/local/seqkit_stats/main.nf'
include { TRIM_AND_CONCAT           } from '../../modules/local/custom/trimandconcat/main.nf'

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
        ch_input
    )
    ch_versions = ch_versions.mix(CREATE_DEMUX_DEPENDENCIES.out.versions)

    // MODULE: Demultiplex
    CUTADAPT (
        ch_raw_data,
        CREATE_DEMUX_DEPENDENCIES.out.fw_index,
        CREATE_DEMUX_DEPENDENCIES.out.rv_index,
        params.ulimit
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    // MODULE: Rename the Cutadapt output files
    RENAME (
        CUTADAPT.out.reads,
        CREATE_DEMUX_DEPENDENCIES.out.sample_rename
    )
    ch_versions = ch_versions.mix(RENAME.out.versions)

    ////ch_original = CUTADAPT.out.reads.collect { tuple -> [tuple[0][0], tuple[1]] }.view() 

    // MODULE: Check stats of reads assigned to samples after demultiplexing
    ASSIGNED_STATS (
        RENAME.out.reads,
        "assigned"
    )
    ch_versions = ch_versions.mix(ASSIGNED_STATS.out.versions)

    // MODULE: Check stats of reads that couldn't be assigned to samples after demultiplexing
    UNKNOWN_STATS (
        RENAME.out.unknown,
        "unknown"
    )
    ch_versions = ch_versions.mix(UNKNOWN_STATS.out.versions)

    // MODULE: Trim leftover primers and concatenate files so that there is only one R1 and one R2 file per sample
    TRIM_AND_CONCAT (
        RENAME.out.reads,
        ch_input
    )
    ch_versions = ch_versions.mix(TRIM_AND_CONCAT.out.versions)

    // MODULE: Check stats after trimming and concatenating files
    FINAL_STATS (
        TRIM_AND_CONCAT.out.reads,
        "final"
    )
    ch_versions = ch_versions.mix(FINAL_STATS.out.versions)

    emit:
    reads           = TRIM_AND_CONCAT.out.reads           // channel: [ val(meta), reads ]
    assigned_stats  = ASSIGNED_STATS.out.stats            // channel: [ val(meta), stats ]
    unknown_stats   = UNKNOWN_STATS.out.stats             // channel: [ val(meta), stats ]
    final_stats     = FINAL_STATS.out.stats               // channel: [ val(meta), stats ]
    raw_stats       = RAW_STATS.out.stats                 // channel: [ val(meta), stats ]
    versions        = ch_versions                         // channel: [ versions.yml ]
}