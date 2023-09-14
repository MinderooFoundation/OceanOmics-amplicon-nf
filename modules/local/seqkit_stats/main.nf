process SEQKIT_STATS {
    tag "$meta"
    label 'process_medium'
    container 'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0'

    input:
    tuple val(meta), path(reads)
    val prefix

    output:
    tuple val(meta), path("*.txt")  , emit: stats
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    seqkit stats -j ${task.cpus} -b ${reads} -a > ${meta}_${prefix}_seqkit_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
