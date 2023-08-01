//
// Create a new channel after demultiplexing in the correct nf-core format
//

include { SAMPLESHEET_CHECK    } from '../../modules/local/samplesheet_check'
include { POSTDEMUX_SAMPSHEET  } from '../../modules/local/custom/postdemuxsampsheet/main.nf'

workflow POSTDEMUX_WORKFLOW {
    take:
    demux_reads 
    samplesheet // file: /path/to/samplesheet.csv
    raw_data

    main:
    ch_versions = Channel.empty()

    // MODULE: Create a new samplesheet with demultiplexed fastq files
    POSTDEMUX_SAMPSHEET (
        samplesheet,
        demux_reads,
        raw_data,
        params.obi3_demux
    )
    ch_versions = ch_versions.mix(POSTDEMUX_SAMPSHEET.out.versions)

    // MODULE: This is needed to create a reads channel that's compatible with nf-core modules
    SAMPLESHEET_CHECK ( POSTDEMUX_SAMPSHEET.out.samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_new_fastq_channel(it) }
        .set { reads }
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    emit:
    reads                                                     // channel: [ val(meta), [ reads ] ]
    missing_samples = POSTDEMUX_SAMPSHEET.out.missing_samples // channel: [ missing_samples.txt ]
    versions = ch_versions                                    // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_new_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    def fastq_meta = []
    meta.id        = row.sample

    if (row.fastq_2.isEmpty()) {
        meta.single_end = true
    } else {
        meta.single_end = false
    }

    // add path(s) of the fastq file(s) to the meta map
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]

    } else {
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }

    return fastq_meta
}
