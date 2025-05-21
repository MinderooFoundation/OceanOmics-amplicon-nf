process CONCATFILE_INFO {
    tag "$prefix"
    label 'process_low'
    container 'ubuntu:20.04'

    input:
    tuple val(prefix), path(csvs)

    output:
    path("*_concat.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo samp_name,lib_id,filename,filename2,checksum_filename,checksum_filename2 > ${prefix}_concat.csv

    cat ${csvs} >> ${prefix}_concat.csv
    """
}
