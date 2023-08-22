process OCOMNBC {
    tag "$prefix"
    label 'process_high'
    container 'adbennett/ocom_nbc:v1.4'

    input:
    tuple val(prefix), path(fasta)

    output:
    tuple val(prefix), path("*.tsv")         , emit: nbc_output
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    """
    naiveBayesClassify.py -i ${fasta} -o ${prefix}_ocom_nbc_output.tsv -t ${task.process} -v
    """
}
