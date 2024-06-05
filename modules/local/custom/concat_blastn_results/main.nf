process CONCAT_BLASTN_RESULTS {
    tag ""
    label 'process_medium'
    container 'quay.io/biocontainers/csvtk:0.26.0--h9ee0642_0'

    input:
    tuple val(prefix), path(files)

    output:
    tuple val(prefix), path("*.txt"), emit: txt
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat $files > ${prefix}_blastn_results.txt

    if [ ! -s "${prefix}_blastn_results.txt" ]; then
        echo "Error: ${prefix} blast results are empty, try a different blast database or try lowering the '--perc_identity' or '--qcov' values."
        exit 1
    fi

    asv_count=\$(cut -f1 "${prefix}_blastn_results.txt" | sort | uniq | wc -l)
    if [ "\$asv_count" -lt 2 ]; then
        echo "Error: blast results have less than 2 ${prefix}s left, try a different blast database or try lowering the '--perc_identity' or '--qcov' values."
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: 4.0
        csvtk: 0.26.0
    END_VERSIONS
    """
}
