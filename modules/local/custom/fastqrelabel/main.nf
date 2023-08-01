process FASTQRELABEL {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container 'ubuntu:20.04'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*relabeled.fq.gz"), emit: reads
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    i=1
    for file in ${reads}
    do
        #zcat "\$file" | sed "s/^@.*/@\${prefix}/" | gzip > "\${prefix}_R\${i}_relabeled.fq.gz"
        zcat "\$file" | awk 'NR % 4 == 1 {print "@'${prefix}'"} NR % 4 != 1 {print}' | gzip > "${prefix}_R\${i}_relabeled.fq.gz"
        i=\$((i+1))                                
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: 1.3.4
    END_VERSIONS
    """
}