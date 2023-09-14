process OBITOOLS3_NGSFILTER {
    tag "$prefix"
    label 'process_medium'
    container 'adbennett/obitools3:v3.0.1b22'

    input:
    tuple val(prefix), path(reads)
    path ngsfile

    output:
    tuple val(prefix), path("*_identified.fq.gz")   , emit: reads
    tuple val(prefix), path("*_unidentified.fq.gz") , emit: unidentified
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    obi import --ngsfilter-input ${ngsfile} view/ngsfile

    obi import --fastq-input --quality-sanger ${reads} view/reads

    obi ngsfilter -t view/ngsfile -u view/unidentified_reads view/reads view/identified_reads

    obi export view/identified_reads --fastq-output -o ${prefix}_identified.fq
    obi export view/unidentified_reads --fastq-output -o ${prefix}_unidentified.fq
    gzip ${prefix}_identified.fq
    gzip ${prefix}_unidentified.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        obitools3: 3.0.1b22
    END_VERSIONS
    """
}
