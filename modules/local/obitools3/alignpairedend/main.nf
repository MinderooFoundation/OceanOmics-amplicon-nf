process OBITOOLS3_ALIGNPAIREDEND {
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
    """
    files=($reads)
    if [ \${#files[@]} -gt 1 ]; then
        i=1
        for file in "\${files[@]}"; do
            obi import --fastq-input --quality-sanger \$file view/read\$i
            i=\$((i + 1))
        done

        obi alignpairedend -R view/read2 view/read1 view/aligned

        obi export view/aligned --fastq-output -o ${prefix}_aligned.fq
        gzip ${prefix}_aligned.fq
    else
        mv \${files[0]} ${prefix}_aligned.fq.gz
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        obitools3: 3.0.1b22
    END_VERSIONS
    """
}
