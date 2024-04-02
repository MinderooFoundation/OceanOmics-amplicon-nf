process CUTADAPT {
    tag "$prefix"
    label 'process_medium'
    container 'quay.io/biocontainers/cutadapt:4.7--py310h4b81fae_1'

    input:
    tuple val(prefix), path(reads)
    path fw_index
    path rv_index
    val ulimit

    output:
    tuple val(prefix), path("*.fq.gz"), emit: reads
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    #!/bin/bash

    # Too avoid too many open files error:
    ulimit -S -n ${ulimit}

    files=($reads)
    if [ "\${#files[@]}" -eq 1 ]; then
        cutadapt -j ${task.cpus} \
            ${args1} \
            ${reads}
    else
        cutadapt -j ${task.cpus} \
            ${args2} \
            ${reads}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
