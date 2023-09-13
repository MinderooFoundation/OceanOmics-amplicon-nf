process SPLIT_FASTQ {
    tag "$samplesheet"
    label 'process_medium'
    container 'quay.io/biocontainers/csvtk:0.26.0--h9ee0642_0'

    input:
    tuple val(prefix), path(reads)
    path samplesheet

    output:
    tuple val(prefix), path("*.fq.gz"), emit: reads
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    zcat ${reads} > tmp.fq

    cat ${samplesheet} | csvtk cut -f sample | {
        skip_first=true  # Flag variable to skip the first iteration

        while IFS= read -r sample; do
            if \$skip_first; then
                skip_first=false
                continue
            fi

            cat tmp.fq | grep -A 3 "sample=\${sample}" | grep -v "^--\$" > \${sample}.R1.fq &
        done
        wait
    }

    rm tmp.fq
    gzip *.R1.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: 0.26.0
    END_VERSIONS
    """
}
