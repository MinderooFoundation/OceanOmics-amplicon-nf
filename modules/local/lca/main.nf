process LCA {
    tag "$prefix"
    label 'process_medium'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(prefix), path(table), path(blast_results)

    output:
    path "*intermediate.tab"                  , emit: intermediate
    tuple val(prefix), path("*lca_output.tab"), emit: lca_output
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    """
    runAssign_collapsedTaxonomy.py \\
        $table \\
        $blast_results \\
        $args \\
        ${prefix}_lca_output.tab
    
    mv interMediate_res.tab ${prefix}_intermediate.tab

    asv_count=\$(tail -n +2 "${prefix}_lca_output.tab" | wc -l)
    if [ "\$asv_count" -lt 2 ]; then
        echo "Error: LCA results have less than 2 ${prefix}s left, try changing the LCA values '--lca_qcov', '--lca_pid', '--lca_diff'."
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
