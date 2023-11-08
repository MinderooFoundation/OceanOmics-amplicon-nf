process SEQTK_TRIM {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/seqtk:1.4--he4a0461_1'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("*.gz"), emit: reads
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for i in $read
    do
        zcat \$i > tmp.fq

        seqtk trimfq $args \$i > tmp.fq

        gzip tmp.fq
        mv tmp.fq.gz trimmed_\${i}
        touch \$i
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: 1.4
    END_VERSIONS
    """
}
