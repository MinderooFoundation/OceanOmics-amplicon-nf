process LULU {
    tag "$prefix"
    label 'process_medium'
    container 'mahsamousavi/lulu:2019'

    input:
    tuple val(prefix), path(table)
    path match_list

    output:
    tuple val(prefix), path("*curated_table.tab"), emit: curated_table
    path "*lulu_map.tab"                         , emit: lulu_map
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "\"${prefix}\""
    def table = "\"${table}\""
    def match_list = "\"${match_list}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(lulu))
    # This is a lulu script for post-clustering of Zotu table created by unoise3

    otutab <- read.table($table, header = TRUE, check.names = FALSE, sep = "\\t", as.is=TRUE, row.names = 1)
    matchlist <- read.table($match_list, header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

    sequence <- data.frame(names = rownames(otutab), sequence = otutab[, c(paste0(toupper($prefix), "_sequence"))])
    otutab[, c(paste0(toupper($prefix), "_sequence"))] <- list(NULL)

    # Curation step with lulu
    curated_result <- lulu(otutab, matchlist, $args) # This runs the default parameter of lulu (i.e. minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)

    # Format table for LCA iput
    curated_table <- curated_result\$curated_table
    curated_table[, paste0(toupper($prefix), "_sequence")] <- sequence\$sequence[match(rownames(curated_table), sequence\$names)]

    write.table(curated_table, paste0($prefix, "_curated_table.tab"), sep="\\t", quote = FALSE, row.names = FALSE)  # write curated result
    write.table(curated_result\$otu_map, paste0($prefix, "_lulu_map.tab"), sep="\\t")            # write the map info

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    lulu: ", packageVersion("lulu"))), "versions.yml")
    """
}
