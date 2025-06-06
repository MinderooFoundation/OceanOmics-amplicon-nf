process FLAG_OTUS_OUTSIDERANGE {
    tag "$prefix"
    label 'process_medium'

    container 'adbennett/phyloseq_and_tree:v2'

    input:
    tuple val(prefix), path(phyloseq)
    val(min_length)
    val(max_length)

    output:
    tuple val(prefix), path("*_flagged_phyloseq.rds"), emit: phyloseq_object

    when:
    task.ext.when == null || task.ext.when

    script:
    def phyloseq   = "\"${phyloseq}\""
    def prefix     = "\"${prefix}\""
    def min_length = "\"${min_length}\""
    def max_length = "\"${max_length}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(phyloseq))

    phyloseq <- readRDS(${phyloseq})
    TAX  <- phyloseq@tax_table@.Data
    OTU  <- phyloseq@otu_table
    SAM  <- phyloseq@sam_data
    TREE <- phyloseq@phy_tree

    if ("ASV_sequence" %in% colnames(TAX)) {
        TAX = cbind(TAX, unusual_size = nchar(TAX[, "ASV_sequence"]) >= ${min_length} | nchar(TAX[, "ASV_sequence"]) <= ${max_length})
    }
    if ("ZOTU_sequence" %in% colnames(TAX)) {
        TAX = cbind(TAX, unusual_size = nchar(TAX[, "ZOTU_sequence"]) >= ${min_length} | nchar(TAX[, "ZOTU_sequence"]) <= ${max_length})
    }

    out_phyloseq <- phyloseq(tax_table(TAX), OTU, SAM, TREE)
    saveRDS(out_phyloseq, paste0(${prefix}, "_flagged_phyloseq.rds"))
    """
}
