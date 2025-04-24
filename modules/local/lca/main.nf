process LCA {
    tag "$prefix"
    label 'process_medium'
    label 'error_retry'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(prefix), path(table), path(blast_results), path(fasta)
    path db

    output:
    path "*intermediate.tab"                  , emit: intermediate
    tuple val(prefix), path("*lca_output.tab"), emit: lca_output
    tuple val(prefix), path("*taxaRaw.tsv")   , emit: taxa_raw
    tuple val(prefix), path("*taxaFinal.tsv") , emit: taxa_final
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    """
    # Create an associative array from the FASTA file
    awk '
    BEGIN { FS = "\\n"; RS = ">" }
    NR > 1 {
        split(\$1, header, " ")
        id = header[1]
        gsub(/\\n/, "", \$2)
        seq = \$2
        fasta[id] = seq
    }
    END {
        for (i in fasta) {
            print i "\\t" fasta[i]
        }
    }' *.fa > id_to_seq.tsv

    #  Join the TSV file with the FASTA sequences using the id (column 1)
    # Make sure both files are sorted by the ID column
    sort -k1,1 *_blastn_results.txt > sorted_BLAST.tsv
    sort -k1,1 id_to_seq.tsv > sorted_id_to_seq.tsv

    # Step 3: Join them together
    awk -F'\\t' 'FNR==NR { seq[\$1] = \$2; next } { s = (\$1 in seq ? seq[\$1] : "NA"); print \$0 "\\t" s }' sorted_id_to_seq.tsv sorted_BLAST.tsv > BLAST_with_seq.tsv

    # Cleanup
    rm sorted_BLAST.tsv sorted_id_to_seq.tsv id_to_seq.tsv

    DB=`ls *.ndb | sed 's/\\.ndb\$//'`

    runAssign_collapsedTaxonomy.py \\
        $table \\
        BLAST_with_seq.tsv \\
        $args \\
        ${prefix}_lca_output.tab \\
        \$DB

    mv interMediate_res.tab ${prefix}_intermediate.tab
    cat taxaRaw.tsv | sed 's/"//g' > ${prefix}_taxaRaw.tsv
    cat taxaFinal.tsv | sed 's/"//g' > ${prefix}_taxaFinal.tsv
    rm taxaRaw.tsv
    rm taxaFinal.tsv

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
