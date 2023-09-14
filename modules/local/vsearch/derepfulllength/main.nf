process VSEARCH_DEREPFULLLENGTH {
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
    vsearch --derep_fulllength ${reads} ${args} --output Unq.fa
    gzip Unq.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
