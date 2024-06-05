process GETLULUFASTA {
    label 'process_low'

    container 'ubuntu:20.04'

    input:
    tuple val(prefix), path(fasta)
    tuple val(prefix), path(curated_table)

    output:
    tuple val(prefix), path("*.fa"), emit: fasta
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    touch curated_$fasta
    cut -f1 $curated_table > ${prefix}_names
    tail -n +2 "${prefix}_names" | while IFS= read -r line
    do
        asv_name=\$(echo \$line | cut -f 1)
        cat $fasta | grep \${line}\$ -A 1 >> curated_$fasta
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
