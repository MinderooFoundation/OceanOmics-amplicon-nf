process LCA_WITH_FISHBASE {
    tag "$prefix"
    label 'process_medium'
    //label 'error_retry'
    //container 'quay.io/biocontainers/python:3.9--1'
    //container 'quay.io/biocontainers/pandas:2.2.1'
    container 'quay.io/microbiome-informatics/pandas-pyarrow:pd2.2.1_pya15.0.0'

    input:
    tuple val(prefix), path(table), path(blast_results)
    path(worms)

    output:
    tuple val(prefix), path("*_lca_with_fishbase_output.tsv"), emit: lca_output
    tuple val(prefix), path("*taxa_raw.tsv")                 , emit: taxa_raw
    tuple val(prefix), path("*taxa_final.tsv")               , emit: taxa_final
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    """
    calculateLCAWithFishbase_Claude_FAIReCompatible.py \\
        -f ${blast_results} \\
        -o ${prefix}_lca_with_fishbase_output.tsv \\
        --asv_table ${table} \\
        --raw_output taxaRaw.tsv \\
        --final_output taxaFinal.tsv \\
        --worms_file ${worms} \\
        $args

    cat taxaRaw.tsv | sed 's/"//g' > ${prefix}_taxa_raw.tsv
    cat taxaFinal.tsv | sed 's/"//g' > ${prefix}_taxa_final.tsv
    rm taxaRaw.tsv
    rm taxaFinal.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
