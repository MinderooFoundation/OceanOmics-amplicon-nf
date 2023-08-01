process OBITOOLS3_GREP {
    tag "$prefix"
    label 'process_medium'
    container 'adbennett/obitools3:v3.0.1b22'

    input:
    tuple val(prefix), path(reads)

    output:
    tuple val(prefix), path("*.fq.gz"), emit: reads
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    obi import --fastq-input --quality-sanger $reads view/reads
    
    obi grep $args view/reads view/filtered_reads

    obi export view/filtered_reads --fastq-output -o ${prefix}.R1.fq
    gzip ${prefix}.R1.fq
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        obitools3: 3.0.1b22
    END_VERSIONS
    """
}
