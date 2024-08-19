process NESTER_FILTER {
    tag "$prefix"
    label 'process_low'
    container 'adbennett/phyloseq_and_tree:v2'

    input:
    tuple val(prefix), path(phyloseq_object), path(final_taxa)

    output:
    path("*_final_taxa_filtered.tsv"), emit: final_taxa
    path("*_nester_stats.txt")       , emit: stats
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

    phyloseq   <- readRDS($phyloseq)
    final_taxa <- read.table($final_taxa, sep="\t")

    TREE         <- phyloseq@phy_tree
    OTU          <- data.frame(phyloseq@otu_table)
    TAX          <- data.frame(phyloseq@tax_table)
    SAM          <- data.frame(phyloseq@sam_data)
    controls     <- rownames(SAM[SAM\$sample_type == "negative_control", ])
    #controls    <- rownames(SAM[SAM\$sample_type == "negative_control" | SAM\$sample_type == "extraction_blank" | SAM\$sample_type == "DI_control", ])
    stats_vector <- c(paste0("Samples identified as controls: ", controls))

    for (asv in rownames(OTU)) {
        control_read_count <- sum(OTU[asv, colnames(OTU) %in% controls])
        sample_read_count  <- sum(OTU[asv, ! colnames(OTU) %in% controls])
        total_read_count   <- control_read_count + sample_read_count
        control_percent    <- (control_read_count/sample_read_count)*100

        stats_vector <- c(stats_vector, "")
        stats_vector <- c(stats_vector, paste0("ASV: ", asv))
        stats_vector <- c(stats_vector, paste0("read count in controls - ", control_read_count))
        stats_vector <- c(stats_vector, paste0("read count in samples - ", sample_read_count))
        stats_vector <- c(stats_vector, paste0("read count total - ", total_read_count))
        stats_vector <- c(stats_vector, paste0("percentage of reads in controls - ", control_percent, "%"))

        # Drop ASVs in more than 0.5% of all reads come from controls
        if (control_percent > 0.5) {
            stats_vector <- c(stats_vector, paste0("percentage greater than 0.5%, dropping ", asv))
            final_taxa   <- final_taxa[final_taxa\$ASV != asv, ]
            OTU          <- OTU[! rownames(OTU) %in% asv, ]
            TAX          <- TAX[! rownames(TAX) %in% asv, ]

        } else if (control_percent > 0) {
            stats_vector  <- c(stats_vector, paste0("filtering ", asv))
            sample_string <- ""
            before_string <- ""
            after_string  <- ""

            for (sample in colnames(OTU)[!colnames(OTU) %in% controls]) {
                curr_count           <- OTU[asv, sample]
                potential_filt_count <- round(curr_count * control_percent)
                filt_count           <- control_read_count

                if (potential_filt_count > control_read_count) {
                    filt_count <- potential_filt_count
                }

                OTU[asv, sample] <- curr_count - filt_count
                if (OTU[asv, sample] < 0) {
                    OTU[asv, sample] <- 0
                }

                sample_string_len <- nchar(sample)
                before_string_len <- nchar(curr_count)
                after_string_len  <- nchar(OTU[asv, sample])
                max_string_len    <- max(c(sample_string_len, before_string_len, after_string_len)) + 1
                sample_string     <- paste0(sample_string, sample, strrep(" ", max_string_len - sample_string_len))
                before_string     <- paste0(before_string, curr_count, strrep(" ", max_string_len - before_string_len))
                after_string      <- paste0(after_string, OTU[asv, sample], strrep(" ", max_string_len - after_string_len))
            }

            stats_vector <- c(stats_vector, "before: ")
            stats_vector <- c(stats_vector, sample_string)
            stats_vector <- c(stats_vector, before_string)
            stats_vector <- c(stats_vector, "after: ")
            stats_vector <- c(stats_vector, sample_string)
            stats_vector <- c(stats_vector, after_string)

        } else {
            stats_vector <- c(stats_vector, paste0("No filtering needed for ", asv))
        }
    }

    stats_conn <- file(paste0($prefix, "_nester_stats.txt"))
    writeLines(stats_vector, stats_conn)
    close(stats_conn)

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    phyloseq: ", packageVersion("phyloseq"))), "versions.yml")
    """
}
