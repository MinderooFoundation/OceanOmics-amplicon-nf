process RENAME_OTURAW {
    tag "$prefix"
    label 'process_low'

    container 'ubuntu:20.04'

    input:
    tuple val(prefix), path(otu)

    output:
    tuple val(prefix), path("raw_*"), emit: otu_raw
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $otu | sed 's/_T1//g' | sed 's/#ID\\t//g' | sed 's/"//g' > raw_${otu}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 1.3.4
    END_VERSIONS
    """
}
