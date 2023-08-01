process FASTQCONCAT {
    tag "$reads"
    label 'process_medium'

    container 'ubuntu:20.04'

    input:
    path reads

    output:
    path "*.fq.gz"     , emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    zcat ${reads} | gzip > "concat_QF_Dmux.fq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
