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

    otutab <- read.table($table, header = TRUE, check.names = FALSE, sep = "\t", as.is=TRUE, row.names = 1)
    matchlist <- read.table($match_list, header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

    otutab[, c(paste0(toupper($prefix), "_sequence"))] <- list(NULL)

    # Curation step with lulu
    curated_result <- lulu(otutab, matchlist, $args) # This runs the default parameter of lulu (i.e. minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)
    # there are more parameters to play with in lulu command. Check lulu paper to understand how they will affect your results 

    curated_result\$curated_table # Curated OTU table
    curated_result\$curated_count # Number of OTUs retained
    curated_result\$curated_otus # IDs of curated OTUs
    curated_result\$discarded_count # OTUs discarded
    curated_result\$otu_map # total - total read count, spread - the number of samples the OTU is present in
                       # parent_id - ID of OTU with which this OTU was merged (or self)
                       # curated - ("parent" or "merged"), was this OTU kept as a valid OTU (parent) or merged with another
                       # rank - The rank of the OTU in terms of decreasing spread and read count

    curated_result\$original_table # Original OTU table

    # Format table for LCA iput
    curated_table <- curated_result\$curated_table
    curated_table\$'#ID' <- rownames(curated_table)
    curated_table <- curated_table[, c(ncol(curated_table), 1:(ncol(curated_table)-1))]

    write.table(curated_table, paste0($prefix, "_curated_table.tab"), sep="\t", quote = FALSE, row.names = FALSE)  # write curated result 
    write.table(curated_result\$otu_map, paste0($prefix, "_lulu_map.tab"), sep="\t")            # write the map info


    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    lulu: ", packageVersion("lulu"))), "versions.yml")
    """
}