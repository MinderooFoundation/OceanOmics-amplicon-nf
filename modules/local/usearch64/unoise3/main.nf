process USEARCH64_UNOISE3 {
    tag "$reads"
    label 'process_medium'
    container 'sunqiangkun/usearch:v1'

    input:
    val usearch64
    val prefix
    path reads 

    output:
    tuple val(prefix), path("*.fa") , emit: zotu_fasta
    tuple val(prefix), path("*.txt"), emit: tabbed_out
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    zcat ${reads} > input.fa
    ${usearch64} -unoise3 input.fa -zotus tmp.fa -tabbedout Unq_unoise3.txt ${args}
    rm input.fa

    # This makes sure each read is only two lines
    awk '/^>/ {if (seq) print seq; printf("%s\\n",\$0); seq=""; next} {seq = seq \$0} END {if (seq) print seq}' tmp.fa > zotus.fa
    rm tmp.fa

    if [ ! -s "zotus.fa" ]; then
        echo "Error: 0 zotus after unoise3, at least 2 are needed. Try changing usearch option values ('--min_size', '--min_quality', or '--min_align_len')."
        exit 1
    fi

    zotus_count=\$(grep -c '^>' "zotus.fa")
    if (( zotus_count < 2 )); then
        echo "Error: \${zotus_count} zotu after unoise3, at least 2 are needed. Try changing usearch option values ('--min_size', '--min_quality', or '--min_align_len')."
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch: \$(usearch --version 2>&1 | head -n 1 | sed 's/usearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}