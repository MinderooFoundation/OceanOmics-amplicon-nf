process ADAPTERREMOVAL {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/adapterremoval:2.3.2--hb7ba0dd_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*fq.gz"), emit: reads
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        AdapterRemoval \\
            --threads ${task.cpus} \\
            --file1 ${reads} \\
            --basename ${prefix} \\
            $args

        # Get rid of the M_ and MT_ AdapterRemoval added to the sample names
        cat ${prefix}.truncated | sed 's/^@M_/@/; s/^@MT_/@/' > ${prefix}_QF.fq
        gzip ${prefix}_QF.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            adapterremoval: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    } else {
        """
        AdapterRemoval \\
            --threads ${task.cpus} \\
            --file1 ${reads[0]} \\
            --file2 ${reads[1]} \\
            --basename ${prefix} \\
            $args2

        # Get rid of the M_ and MT_ AdapterRemoval added to the sample names
        cat ${prefix}.collapsed | sed 's/^@M_/@/; s/^@MT_/@/' > ${prefix}_QF.fq
        gzip ${prefix}_QF.fq

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            adapterremoval: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    }
}
