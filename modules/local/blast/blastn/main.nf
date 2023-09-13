process BLAST_BLASTN {
    tag "$prefix"
    label 'process_high'

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
    DB=`ls *.ndb | sed 's/\\.ndb\$//'`
    blastn \\
        -num_threads $task.cpus \\
        -db \$DB \\
        -query $fasta \\
        $args \\
        -out ${prefix}_blastn_results.txt

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
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
