process VSEARCH_CLUSTERUNOISE {
    tag "$reads"
    label 'process_medium'
    container 'quay.io/biocontainers/vsearch:2.21.1--h95f258a_0'

    input:
    path reads

    output:
    path "*.fa.gz"     , emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vsearch --cluster_unoise ${reads} ${args} centroids.fa

    if [ ! -s "centroids.fa" ]; then
        echo "Error: 0 zotus after cluster_unoise, at least 2 are needed. Try changing vsearch option values ('--min_size', '--min_quality', or '--min_align_len')."
        exit 1
    fi

    zotus_count=\$(grep -c '^>' "centroids.fa")
    if (( zotus_count < 2 )); then
        echo "Error: \${zotus_count} zotu after cluster_unoise, at least 2 are needed. Try changing vsearch option values ('--min_size', '--min_quality', or '--min_align_len')."
        exit 1
    fi

    gzip centroids.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
