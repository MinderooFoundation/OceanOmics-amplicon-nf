process CURATE_BLASTN_RESULTS {
    tag ""
    label 'process_medium'
    container 'ubuntu:20.04'

    input:
    tuple val(prefix) , path(fasta)
    tuple val(prefix2), path(blastn_results)

    output:
    tuple val(prefix), path("curated_*"), emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    grep "^>" "${fasta}" | sed 's/^>//' > asv_names.txt
    grep -wf asv_names.txt "${blastn_results}" > "curated_${blastn_results}"
    """
}
