process FASTQTOFASTA {
    tag "$reads"
    label 'process_low'

    container 'ubuntu:20.04'

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
    zcat ${reads} | sed -n '1~4s/^@/>/p;2~4p' | gzip > concat.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
