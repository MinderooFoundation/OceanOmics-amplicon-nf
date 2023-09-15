/*
 * Demultiplex with Obitool3
 */

include { CREATE_NGSFILE           } from '../../modules/local/custom/createngsfile/main.nf'
include { OBITOOLS3_ALIGNPAIREDEND } from '../../modules/local/obitools3/alignpairedend/main.nf'
include { OBITOOLS3_NGSFILTER      } from '../../modules/local/obitools3/ngsfilter/main.nf'
include { OBITOOLS3_GREP           } from '../../modules/local/obitools3/grep/main.nf'
include { SPLIT_FASTQ              } from '../../modules/local/custom/splitfastq/main.nf'
include { SEQKIT_STATS as \
            ASSIGNED_STATS;
            SEQKIT_STATS as \
            UNKNOWN_STATS;
            SEQKIT_STATS as \
            RAW_STATS;
            SEQKIT_STATS as \
            FINAL_STATS            } from '../../modules/local/seqkit_stats/main.nf'

workflow OBITOOLS3_WORKFLOW {
    take:
    ch_input    // channel: [ samplesheet.csv ]
    ch_raw_data // channel: [ prefix, [ reads ] ]

    main:
    ch_versions = Channel.empty()

    // MODULE: Check stats of raw data
    RAW_STATS (
        ch_raw_data,
        "raw"
    )
    ch_versions = ch_versions.mix(RAW_STATS.out.versions)

    // MODULE: Align if paired end
    OBITOOLS3_ALIGNPAIREDEND (
        ch_raw_data
    )
    ch_versions = ch_versions.mix(OBITOOLS3_ALIGNPAIREDEND.out.versions)

    // MODULE: Create NGS file for demultiplexing
    CREATE_NGSFILE (
        ch_input
    )
    ch_versions = ch_versions.mix(CREATE_NGSFILE.out.versions)

    // MODULE: Demultiplex
    OBITOOLS3_NGSFILTER (
        OBITOOLS3_ALIGNPAIREDEND.out.reads,
        CREATE_NGSFILE.out.ngsfile
    )
    ch_versions = ch_versions.mix(OBITOOLS3_NGSFILTER.out.versions)

    // MODULE: Check stats of reads assigned to samples after demultiplexing
    ASSIGNED_STATS (
        OBITOOLS3_NGSFILTER.out.reads,
        "assigned"
    )
    ch_versions = ch_versions.mix(ASSIGNED_STATS.out.versions)

    // MODULE: Check stats of reads that couldn't be assigned to samples after demultiplexing
    UNKNOWN_STATS (
        OBITOOLS3_NGSFILTER.out.unidentified,
        "unknown"
    )
    ch_versions = ch_versions.mix(UNKNOWN_STATS.out.versions)

    // MODULE: Length filter
    OBITOOLS3_GREP (
        OBITOOLS3_NGSFILTER.out.reads
    )
    ch_versions = ch_versions.mix(OBITOOLS3_GREP.out.versions)

    // MODULE: Check stats after filtering files
    FINAL_STATS (
        OBITOOLS3_GREP.out.reads,
        "final"
    )
    ch_versions = ch_versions.mix(FINAL_STATS.out.versions)

    // MODULE: Split samples into seperate fq files
    SPLIT_FASTQ (
        OBITOOLS3_GREP.out.reads,
        ch_input
    )
    ch_versions = ch_versions.mix(SPLIT_FASTQ.out.versions)

    emit:
    reads           = SPLIT_FASTQ.out.reads               // channel: [ val(meta), reads ]
    assigned_stats  = ASSIGNED_STATS.out.stats            // channel: [ val(meta), stats ]
    unknown_stats   = UNKNOWN_STATS.out.stats             // channel: [ val(meta), stats ]
    final_stats     = FINAL_STATS.out.stats               // channel: [ val(meta), stats ]
    raw_stats       = RAW_STATS.out.stats                 // channel: [ val(meta), stats ]
    versions        = ch_versions                         // channel: [ versions.yml ]
}
