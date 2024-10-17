process FASTP {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/fastp:0.23.4--h125f33a_5'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.gz"), emit: reads
    path "*.json"                , emit: json
    path "*.html"                , emit: html
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        fastp \
            --in1 $reads \
            --out1 ${prefix}.fq.gz \
            --json ${prefix}.json \
            --html ${prefix}.html \
            --thread ${task.cpus} \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version | sed "s/fastp //g")
        END_VERSIONS
        """
    } else {
        """
        fastp \
            --in1 ${reads[0]} \
            --in2 ${reads[1]} \
            --out1 ${prefix}.R1.fq.gz \
            --out2 ${prefix}.R2.fq.gz \
            --json ${prefix}.json \
            --html ${prefix}.html \
            --thread ${task.cpus} \
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version | sed "s/fastp //g")
        END_VERSIONS
        """
    }
}
