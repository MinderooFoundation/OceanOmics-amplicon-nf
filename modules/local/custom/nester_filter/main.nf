process NESTER_FILTER {
    tag "$prefix"
    label 'process_low'
    container 'adbennett/phyloseq_and_tree:v2'

    input:
    tuple val(prefix), path(phyloseq_object), path(final_taxa)

    output:
    path("*_final_taxa_filtered*.tsv"), emit: final_taxa
    path("*_nester_stats*.txt")       , emit: stats
    path("*phyloseq_filtered*.rds")   , emit: phyloseq_object
    path "versions.yml"               , emit: versions

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
    suppressPackageStartupMessages(library(DECIPHER))
    suppressPackageStartupMessages(library(phangorn))
    suppressPackageStartupMessages(library(Biostrings))
    suppressPackageStartupMessages(library(stringr))

    phyloseq   <- readRDS($phyloseq_object)
    final_taxa <- read.table($final_taxa, sep="\t")

    OTU          <- data.frame(phyloseq@otu_table, check.names = FALSE)
    TAX          <- data.frame(phyloseq@tax_table)
    SAM          <- data.frame(phyloseq@sam_data)

    if (! "use_for_filter" %in% colnames(SAM)) {
        SAM\$use_for_filter <- FALSE
    } else {
        SAM\$use_for_filter <- as.logical(SAM\$use_for_filter)
    }

    controls     <- rownames(SAM[SAM\$use_for_filter, ])

    stats_vector <- c("Samples identified as controls: ", controls)
    upper_prefix <- toupper($prefix)
    upper_prefix <- str_split_1(upper_prefix, "_")[1]

    for (asv in rownames(OTU)) {
        control_read_count <- sum(OTU[asv, colnames(OTU) %in% controls])
        sample_read_count  <- sum(OTU[asv, ! colnames(OTU) %in% controls])
        total_read_count   <- control_read_count + sample_read_count
        control_percent    <- (control_read_count/total_read_count)*100

        stats_vector <- c(stats_vector, "")
        stats_vector <- c(stats_vector, paste0("ASV: ", asv))
        stats_vector <- c(stats_vector, paste0("read count in controls - ", control_read_count))
        stats_vector <- c(stats_vector, paste0("read count in samples - ", sample_read_count))
        stats_vector <- c(stats_vector, paste0("read count total - ", total_read_count))
        stats_vector <- c(stats_vector, paste0("percentage of reads in controls - ", control_percent, "%"))

        # Drop ASVs in more than 0.5% of all reads come from controls
        if (control_percent > 100) {
            stats_vector <- c(stats_vector, paste0("percentage greater than 0.5%, dropping ", asv))
            final_taxa   <- final_taxa[final_taxa\$ASV != asv, ]
            OTU          <- OTU[! rownames(OTU) %in% asv, ]
            TAX          <- TAX[! rownames(TAX) %in% asv, ]

        } else if (control_percent > 0) {
            stats_vector  <- c(stats_vector, paste0("filtering ", asv))
            sample_string <- ""
            before_string <- ""
            after_string  <- ""

            for (sample in colnames(OTU)) {
                curr_count           <- OTU[asv, sample]
                potential_filt_count <- round(curr_count * (control_percent/100))
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

    write.table(final_taxa, paste0($prefix, "_final_taxa_filtered.tsv"), sep = "\t")

    #if (length(TAX[[paste0(upper_prefix, "_sequence")]]) >= 3) {
        # Phylogenetic tree code based on code from
        # https://ucdavis-bioinformatics-training.github.io/2021-May-Microbial-Community-Analysis/data_reduction/02-dada2
    #    DNA_set = DNAStringSet(TAX[[paste0(upper_prefix, "_sequence")]])
    #    names(DNA_set) = paste0(TAX[[upper_prefix]])

    #    alignment = AlignSeqs(DNA_set, anchor=NA, processors=\${task.cpus})

    #    phang_align <- phyDat(as(alignment, "matrix"), type="DNA")
    #    dm          <- dist.ml(phang_align)
    #    treeNJ      <- NJ(dm)

    #    fit    <- pml(treeNJ, data=phang_align)
    #    fitGTR <- update(fit, k=4, inv=0.2)
    #}

    TAX\$numberOfUnq_BlastHits <- as.numeric(TAX\$numberOfUnq_BlastHits)
    NEW_OTU                    <- otu_table(OTU, taxa_are_rows = TRUE)
    NEW_TAX                    <- tax_table(TAX)
    rownames(NEW_TAX)          <- rownames(TAX)
    colnames(NEW_TAX)          <- colnames(TAX)
    NEW_SAM                    <- sample_data(SAM)

    #if (length(TAX[[paste0(upper_prefix, "_sequence")]]) >= 3) {
    #    TREE   <- phy_tree(fitGTR\$tree)
    #    physeq <- phyloseq(NEW_OTU, NEW_TAX, NEW_SAM, TREE)
    #} else {
        physeq <- phyloseq(NEW_OTU, NEW_TAX, NEW_SAM)
    #}

    saveRDS(physeq, file = paste0($prefix, "_phyloseq_filtered.rds"))

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    phyloseq: ", packageVersion("phyloseq"))), "versions.yml")
    """
}
