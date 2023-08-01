process USEARCH32_FASTXUNIQUES {
    tag "$reads"
    label 'process_medium'
    container 'sunqiangkun/usearch:v1'

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
    zcat ${reads} > input.fa
    usearch -fastx_uniques input.fa ${args} -fastaout Unq.fa
    gzip Unq.fa
    rm input.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch: \$(usearch --version 2>&1 | head -n 1 | sed 's/usearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}