process BLAST_MAKEBLASTDB {
    tag "$prefix"
    label 'process_medium'

    container'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

    input:
    tuple val(prefix), path(fasta)

    output:
    path '*_db'        , emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    makeblastdb \\
        -in $fasta \\
        -out ${prefix} \\
        $args
    mkdir blast_db
    mv ${prefix}* blast_db
    mv blast_db ${prefix}_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
