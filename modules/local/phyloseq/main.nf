process PHYLOSEQ {
    tag "$prefix"
    label 'process_medium'
    container 'adbennett/phyloseq_readr:v1.0'

    input:
    tuple val(prefix), path(otu_table), path(lca_table), path(nbc_table)
    path metadata
    path filter_table

    output:
    tuple val(prefix), path("*phyloseq.rds")   , emit: phyloseq_object
    tuple val(prefix), path("*_final_taxa.tsv"), emit: final_taxa
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def otu_table    = "\"${otu_table}\""
    def lca_table    = "\"${lca_table}\""
    def nbc_table    = "\"${nbc_table}\""
    def metadata     = "\"${metadata}\""
    def filter_table = "\"${filter_table}\""
    def prefix       = "\"${prefix}\""
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(phyloseq))
    suppressPackageStartupMessages(library(readr))
    #suppressPackageStartupMessages(library(DECIPHER))
    #suppressPackageStartupMessages(library(phangorn))
    #suppressPackageStartupMessages(library(Biostrings))
    suppressPackageStartupMessages(library(stringr))

    #........................................................................
    # Prepare merged data
    #........................................................................

    lca_tab = tryCatch({
        read.table($lca_table, sep="\\t", header=TRUE, check.names = FALSE)
    }, error = function(e) {
        lca_tab <- read_tsv($lca_table)
    })
    otu_tab <- read.table($otu_table, sep="\\t", header=TRUE, check.names = FALSE)
    nbc_tab <- read.table($nbc_table, sep="\\t", header=TRUE)
    meta    <- read.csv($metadata, header=TRUE)

    upper_prefix                  <- toupper($prefix)
    upper_prefix                  <- str_split_1(upper_prefix, "_")[1]
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
    rownames(meta)   <- meta\$samp_name
    meta\$samp_name     <- NULL


    #........................................................................
    # Prepare taxa data
    #........................................................................

    taxa["LCA"] <- ""
    taxa <- taxa[, c(upper_prefix, "domain", "phylum", "class", "order", "family", "genus", "species", "LCA", paste0($prefix, "_length"), "numberOfUnq_BlastHits", "%ID", "Gene", "Genus.prediction", "Genus.score", "Species.prediction", "Species.score", paste0(upper_prefix, "_sequence"))]

    taxa     <- as.data.frame(taxa)

    # Add LCA column
    levels <- c("species", "genus", "family", "order", "class", "phylum", "domain")
    for (row in 1:nrow(taxa)) {
        LCA      <- NA
        level_ID <- 1
        while(LCA == "dropped" | LCA == "" | is.na(LCA)) {
            LCA      <- taxa[row, levels[level_ID]]
            level_ID <- level_ID + 1

            if (level_ID == 9) {
                LCA <- "NA"
            }
        }
        taxa[row, "LCA"] <- LCA
    }

    write.table(taxa, file = paste0($prefix, "_final_taxa.tsv"), sep = "\\t", quote=FALSE)

    rownames(taxa)              <- taxa[[upper_prefix]]
    taxa[[upper_prefix]]        <- NULL
    taxa\$numberOfUnq_BlastHits <- as.numeric(taxa\$numberOfUnq_BlastHits)
    taxa                        <- as.matrix(taxa)


    #........................................................................
    # Prepare otu data
    #........................................................................

    otu           <- as.data.frame(otu)
    rownames(otu) <- otu[[upper_prefix]]
    otu[,c(upper_prefix, 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'queryCoverage', paste0($prefix, "_length"), 'numberOfUnq_BlastHits', '%ID', 'Gene', 'Genus.prediction', 'Genus.score', 'Species.prediction', 'Species.score', paste0(upper_prefix, "_sequence"), 'species_in_LCA', 'sources')] <- list(NULL)

    #........................................................................
    # Create phylogenetic tree
    #........................................................................

    #if (length(seq_tab[[paste0(upper_prefix, "_sequence")]]) >= 3) {
        ## Phylogenetic tree code based on code from
        ## https://ucdavis-bioinformatics-training.github.io/2021-May-Microbial-Community-Analysis/data_reduction/02-dada2
        #DNA_set = DNAStringSet(seq_tab[[paste0(upper_prefix, "_sequence")]])
        #names(DNA_set) = paste0(seq_tab[[upper_prefix]])

        #alignment = AlignSeqs(DNA_set, anchor=NA, processors=\${task.cpus})

        #phang_align <- phyDat(as(alignment, "matrix"), type="DNA")
        #dm          <- dist.ml(phang_align)
        #treeNJ      <- NJ(dm)

        #fit    <- pml(treeNJ, data=phang_align)
        #fitGTR <- update(fit, k=4, inv=0.2)
    #}


    #........................................................................
    # Create phyloseq object
    #........................................................................

    cols_to_drop <- c(
        "sequencing_run", "assay",
        "fw_no", "rv_no",
        "fw_index", "rv_index",
        "fw_primer", "rv_primer",
        "fastq_1", "fastq_2",
        "plate", "well"
    )
    meta_tmp <- meta
    for (col in cols_to_drop) {
        if(col %in% colnames(meta)) {
            meta <- meta[ , !(names(meta) %in% col)]
        }
    }
    # This is to avoid a bug in cases where meta has one or zero columns
    if (! is.data.frame(meta)) {
        meta <- meta_tmp
    }

    OTU    <- otu_table(otu, taxa_are_rows = TRUE)
    TAX    <- tax_table(taxa)
    META   <- sample_data(meta)

    if (length(seq_tab[[paste0(upper_prefix, "_sequence")]]) >= 3) {
        #TREE   <- phy_tree(fitGTR\$tree)
        #physeq <- phyloseq(OTU, TAX, META, TREE)
        physeq <- phyloseq(OTU, TAX, META)
    } else {
        physeq <- phyloseq(OTU, TAX, META)
    }

    saveRDS(physeq, file = paste0($prefix, "_phyloseq.rds"))

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    phyloseq: ", packageVersion("phyloseq"))), "versions.yml")
    """
}
