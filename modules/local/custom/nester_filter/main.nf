process NESTER_FILTER {
    tag "$prefix"
    label 'process_low'
    container 'adbennett/phyloseq_and_tree:v2'

    input:
    tuple val(prefix), path(phyloseq_object), path(final_taxa)

    output:
    path("*_final_taxa_filtered.tsv"), emit: final_taxa
    path("*_nester_stats.tsv")       , emit: stats
    path("*phyloseq_filtered.rds")   , emit: phyloseq_object
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = "\"${prefix}\""
    def phyloseq_object = "\"${phyloseq_object}\""
    def final_taxa      = "\"${final_taxa}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(phyloseq))



    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    phyloseq: ", packageVersion("phyloseq"))), "versions.yml")
    """
}
