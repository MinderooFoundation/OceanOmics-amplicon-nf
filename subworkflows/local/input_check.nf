//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    def fastq_meta = []
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()
    boolean skip_demux = params.skip_demux.toBoolean()

    // add index/primer info of the sample(s) to the meta map
    if (!skip_demux) {
        assert row.fw_index != null : "ERROR: Please check input samplesheet -> 'fw_index' is mandatory if not using '--skip_demux' option!\n${meta.id} is missing 'fw_index'"
        assert row.fw_primer != null : "ERROR: Please check input samplesheet -> 'fw_primer' is mandatory if not using '--skip_demux' option!\n${meta.id} is missing 'fw_primer'"

        meta.fw_index   = row.fw_index
        meta.fw_primer  = row.fw_primer

        if (meta.single_end) {
            meta.rv_index   = null
            meta.rv_primer  = null

        } else {
            assert row.rv_index != null : "ERROR: Please check input samplesheet -> 'rv_index' is mandatory if not using '--skip_demux' option!\n${meta.id} is missing 'rv_index'"
            assert row.rv_primer != null : "ERROR: Please check input samplesheet -> 'rv_primer' is mandatory if not using '--skip_demux' option!\n${meta.id} is missing 'rv_primer'"
            
            meta.rv_index   = row.rv_index
            meta.rv_primer  = row.rv_primer
        }

        return meta

    } else {
        // add path(s) of the fastq file(s) to the meta map
        assert file(row.fastq_1).exists() : "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"

        if (meta.single_end) {
            fastq_meta = [ meta, [ file(row.fastq_1) ] ]

        } else {
            assert file(row.fastq_2).exists() : "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
            
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        }

        return fastq_meta
    }
}
