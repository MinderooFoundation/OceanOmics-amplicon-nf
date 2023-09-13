process CUTADAPT {
    tag "$prefix"
    label 'process_medium'
    container 'quay.io/biocontainers/cutadapt:4.1--py37h8902056_1'

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
    def args = task.ext.args ?: ''
    """
    #!/bin/bash

    # Too avoid too many open files error:
    ulimit -S -n ${ulimit}

    files=($reads)
    if [ "\${#files[@]}" -eq 1 ]; then
        cutadapt -j ${task.cpus} \
            ${args} \
            -g ^file:fw.fa  \
            -o {name}.R1.fq.gz \
            ${reads}
    else
        cutadapt -j ${task.cpus} \
            ${args} \
            -g ^file:fw.fa  \
            -G ^file:rv.fa \
            -o {name1}-{name2}.R1.fq.gz \
            -p {name1}-{name2}.R2.fq.gz \
            ${reads}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
