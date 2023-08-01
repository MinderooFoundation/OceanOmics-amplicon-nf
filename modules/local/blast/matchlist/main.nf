process BLAST_MATCHLIST {
    tag "$prefix"
    label 'process_medium'

    container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

    input:
    tuple val(prefix), path(fasta)
    path db

    output:
    path "*.txt"       , emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    """
    blastn \\
        -num_threads $task.cpus \\
        -db $db/$prefix \\
        -query $fasta \\
        $args \\
        -out ${prefix}_match_list.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}