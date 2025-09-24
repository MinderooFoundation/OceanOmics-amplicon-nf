process BLAST_BLASTN {
    tag "$prefix"
    label 'process_medium'

    container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

    input:
    tuple val(prefix), path(fasta)
    path db

    output:
    tuple val(prefix), path("*.txt"), emit: txt
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    """
    DB=`echo \$(ls *.ndb | sed 's/\\.ndb\$//')`

    blastn \\
        -num_threads $task.cpus \\
        -db "\$DB" \\
        -query $fasta \\
        $args \\
        -out ${prefix}_${fasta}_blastn_results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
