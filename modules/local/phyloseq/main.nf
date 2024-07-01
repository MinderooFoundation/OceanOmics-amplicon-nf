process PHYLOSEQ {
    tag "$prefix"
    label 'process_medium'
    container 'adbennett/phyloseq_and_tree:v2'

    input:
    tuple val(prefix), path(otu_table), path(lca_table), path(nbc_table)
    path metadata
    path filter_table

    output:
    path "*phyloseq.rds"   , emit: phyloseq_object
    path "*_final_taxa.tsv", emit: final_taxa
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def otu_table = "\"${otu_table}\""
    def lca_table = "\"${lca_table}\""
    def nbc_table = "\"${nbc_table}\""
    def metadata = "\"${metadata}\""
    def filter_table = "\"${filter_table}\""
    def prefix = "\"${prefix}\""
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(phyloseq))
    suppressPackageStartupMessages(library(DECIPHER))
    suppressPackageStartupMessages(library(phangorn))
    suppressPackageStartupMessages(library(Biostrings))

    #........................................................................
    # Prepare merged data
    #........................................................................

    lca_tab <- read.table($lca_table, sep="\\t", header=TRUE, check.names = FALSE)
    otu_tab <- read.table($otu_table, sep="\\t", header=TRUE, check.names = FALSE)
    nbc_tab <- read.table($nbc_table, sep="\\t", header=TRUE)
    meta    <- read.table($metadata, sep=",", header=TRUE)

    upper_prefix                  <- toupper($prefix)
    colnames(otu_tab)[1]          <- upper_prefix
    colnames(lca_tab)[8]          <- upper_prefix
    colnames(nbc_tab)[2]          <- upper_prefix
    colnames(lca_tab)             <- sub("_T\\\\d+\$", "", colnames(lca_tab))

    otu_tab <- merge(otu_tab, nbc_tab, by = upper_prefix)

    # Filter if filter table has been provided
    if (file.exists($filter_table)) {
        filt_tab <- read.table($filter_table, sep=",", header=TRUE, na.strings = c("NA", "na", "Na", "nA"))

        for (i in 1:nrow(filt_tab)) {
            level <- filt_tab[i, "level"]
            name  <- filt_tab[i, "name"]

            if (is.na(name)) {
                lca_tab <- subset(lca_tab, !(is.na(lca_tab[[level]])))
            } else {
                lca_tab <- subset(lca_tab, !(lca_tab[[level]] == name) | is.na(lca_tab[[level]]))
            }
        }
    }

    # Add Contam column if control column is in metadata
    if ("control" %in% colnames(meta)) {
        lca_tab\$contam <- "False"

        for (i in seq_along(meta\$sample)) {
            curr_sample <- meta\$sample[i]
            is_control  <- toupper(meta\$control[i])

            if (is_control == "TRUE") {
                lca_tab\$contam[lca_tab[curr_sample] > 0] <- "True"
            }
        }
    }

    otu_nbc    <- otu_tab[, c(upper_prefix, paste0(upper_prefix, "_sequence"), "Gene", "Genus.prediction", "Genus.score", "Species.prediction", "Species.score")]

    # Get the ASVs that are missing blast results
    otu_tab_asvs <- otu_tab[[upper_prefix]]
    lca_tab_asvs <- lca_tab[[upper_prefix]]
    missing_asvs <- setdiff(otu_tab_asvs, lca_tab_asvs)

    df_missing_asvs <- otu_tab[otu_tab[[upper_prefix]] %in% missing_asvs, ]
    if (nrow(df_missing_asvs) > 0) {
        for (i in 1:nrow(df_missing_asvs)) {
            curr_index = nrow(lca_tab) + 1
            lca_tab[curr_index, ] <- NA
            for (y in colnames(df_missing_asvs)) {
                value <- df_missing_asvs[i, y]
                reformated_col <- sub("_T1", "", y)

                if (reformated_col %in% colnames(lca_tab)) {
                    lca_tab[curr_index, reformated_col] <- value
                }
            }
        }
    }

    merged_tab           <- merge(lca_tab, otu_nbc, by = upper_prefix, all.x = TRUE)
    colnames(merged_tab) <- sub(" ", ".", colnames(merged_tab))
    taxa                 <- merged_tab
    otu                  <- merged_tab
    seq_tab              <- merged_tab

    #........................................................................
    # Prepare metadata
    #........................................................................

    meta             <- as.data.frame(meta)
    rownames(meta)   <- meta\$sample
    meta\$sample     <- NULL


    #........................................................................
    # Prepare taxa data
    #........................................................................

    taxa["LCA"] <- ""
    if ("control" %in% colnames(meta)) {
        taxa <- taxa[, c(upper_prefix, "domain", "phylum", "class", "order", "family", "genus", "species", "LCA", "Gene", "Genus.prediction", "Genus.score", "Species.prediction", "Species.score", "contam", paste0(upper_prefix, "_sequence"))]
    } else {
        taxa <- taxa[, c(upper_prefix, "domain", "phylum", "class", "order", "family", "genus", "species", "LCA", "Gene", "Genus.prediction", "Genus.score", "Species.prediction", "Species.score", paste0(upper_prefix, "_sequence"))]
    }

    taxa     <- as.data.frame(taxa)

    # Add LCA column
    levels <- c("species", "genus", "family", "order", "class", "phylum", "domain")
    for (row in 1:nrow(taxa)) {
        LCA      <- NA
        level_ID <- 1
        while(LCA == "dropped" | is.na(LCA)) {
            LCA      <- taxa[row, levels[level_ID]]
            level_ID <- level_ID + 1

            if (level_ID == 9) {
                LCA <- "NA"
            }
        }
        taxa[row, "LCA"] <- LCA
    }

    write.table(taxa, file = paste0($prefix, "_final_taxa.tsv"), sep = "\\t")

    rownames(taxa)       <- taxa[[upper_prefix]]
    taxa[[upper_prefix]] <- NULL
    taxa                 <- as.matrix(taxa)


    #........................................................................
    # Prepare otu data
    #........................................................................

    otu           <- as.data.frame(otu)
    rownames(otu) <- otu[[upper_prefix]]
    otu[,c(upper_prefix, 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'Gene', 'Genus.prediction', 'Genus.score', 'Species.prediction', 'Species.score', 'contam', paste0(upper_prefix, "_sequence"))] <- list(NULL)


    #........................................................................
    # Create phylogenetic tree
    #........................................................................

    if (length(seq_tab[[paste0(upper_prefix, "_sequence")]]) >= 3) {
        # Phylogenetic tree code based on code from
        # https://ucdavis-bioinformatics-training.github.io/2021-May-Microbial-Community-Analysis/data_reduction/02-dada2
        DNA_set = DNAStringSet(seq_tab[[paste0(upper_prefix, "_sequence")]])
        names(DNA_set) = paste0(seq_tab[[upper_prefix]])

        alignment = AlignSeqs(DNA_set, anchor=NA, processors=${task.cpus})

        phang_align <- phyDat(as(alignment, "matrix"), type="DNA")
        dm          <- dist.ml(phang_align)
        treeNJ      <- NJ(dm)

        fit    <- pml(treeNJ, data=phang_align)
        fitGTR <- update(fit, k=4, inv=0.2)
    }


    #........................................................................
    # Create phyloseq object
    #........................................................................

    OTU    <- otu_table(otu, taxa_are_rows = TRUE)
    TAX    <- tax_table(taxa)
    META   <- sample_data(meta)

    if (length(seq_tab[[paste0(upper_prefix, "_sequence")]]) >= 3) {
        TREE   <- phy_tree(fitGTR\$tree)
        physeq <- phyloseq(OTU, TAX, META, TREE)
    } else {
        physeq <- phyloseq(OTU, TAX, META)
    }

    saveRDS(physeq, file = paste0($prefix, "_phyloseq.rds"))

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    phyloseq: ", packageVersion("phyloseq")),paste0("    DECIPHER: ", packageVersion("DECIPHER")),paste0("    phangorn: ", packageVersion("phangorn")),paste0("    Biostrings: ", packageVersion("Biostrings"))), "versions.yml")
    """
}
