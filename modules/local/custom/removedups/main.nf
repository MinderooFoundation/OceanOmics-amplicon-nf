process REMOVE_DUPS {
    tag "$prefix"
    label 'process_medium'

    container 'ubuntu:20.04'

    input:
    tuple val(prefix), path(lca_output)

    output:
    tuple val(prefix), path("*.tsv"), emit: tsv
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -F'\\t' '!seen[\$8]++' ${lca_output} > ${prefix}_lca_deduped.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 1.3.4
    END_VERSIONS
    """
}
