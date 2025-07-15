process ADD_EXPERIMENTRUNMETADATA {
    tag "$prefix"
    label 'process_low'
    //container 'docker.io/pedrofeijao/pandas-openpyxl:v1.0'
    //container 'docker.io/pawelqs/tidyverse_jsonlite_openxlsx:v1'
    //container 'docker.io/adbennett/phyloseq_and_tree:v2'
    container 'docker.io/adbennett/phyloseq_openxlsx_readxl:v1.0'
    //container 'quay.io/biocontainers/bioconductor-phyloseq'

    input:
    //tuple path(phyloseq), path()
    tuple val(prefix), path(phyloseq)
    path(metadata)
    path(samplesheet)
    tuple val(var), path(prefilter_stats)
    val(assay)
    path(inputfile_info)

    output:
    path "*xlsx"       , emit: xlsx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "\"${prefix}\""
    def metadata = "\"${metadata}\""
    def samplesheet = "\"${samplesheet}\""
    def phyloseq = "\"${phyloseq}\""
    def prefilter_stats = "\"${prefilter_stats}\""
    def assay = "\"${assay}\""
    def inputfile_info = "\"${inputfile_info}\""
    """
    #!/usr/bin/env Rscript
    # Get R version
    r_version <- R.version.string

    library(readxl)
    library(openxlsx)
    library(phyloseq)

    getLastRow <- function(wb, sheet) {
        nrows <- nrow(readWorkbook(wb, sheet = sheet))
        if (is.na(nrows)) return(1)
        return(nrows + 1)
    }

    # Load existing workbook
    wb <- loadWorkbook(${metadata})

    phyloseq     <- readRDS(${phyloseq})
    phy_otu      <- data.frame(phyloseq@otu_table)
    phy_tax      <- data.frame(phyloseq@tax_table)
    fileinfo     <- read.csv(${inputfile_info}, header=TRUE, quote="")
    sampsheet    <- read.csv(${samplesheet}, header=TRUE, quote="")
    seqkit_stats <- read.table(${prefilter_stats}, header=TRUE, quote="")
    faire_meta   <- read_excel(${metadata}, sheet="sampleMetadata")
    assay        <- ""

    # Skip non header rows
    continue = TRUE
    if (! "samp_name" %in% colnames(faire_meta)) {
        if (! "samp_name" %in% c(faire_meta[, 1])[[1]]) {
            # Save workbook
            saveWorkbook(wb, paste0(${prefix}, "_faire_metadata.xlsx"), overwrite = TRUE)
            continue = FALSE
        }

        rows_to_skip <- which(faire_meta[, 1] == "samp_name")
        faire_meta   <- read_excel(${metadata}, sheet="sampleMetadata", skip=rows_to_skip)
    }

    if (continue) {
        # If paired-end, remove reverse read rows from stats
        seqkit_stats <- seqkit_stats[seq(1, nrow(seqkit_stats), by = 2), ]

        fileinfo     <- fileinfo[order(fileinfo\$samp_name), ]
        sampsheet    <- sampsheet[order(sampsheet\$samp_name), ]
        seqkit_stats <- seqkit_stats[order(seqkit_stats\$file), ]
        faire_meta   <- faire_meta[faire_meta\$samp_name %in% sampsheet\$samp_name, ]
        new_rows <- data.frame(samp_name = sampsheet[! sampsheet\$samp_name %in% faire_meta\$samp_name, "samp_name"])

        faire_meta = tryCatch({
            new_rows[setdiff(names(faire_meta), names(new_rows))] <- NA
            rbind(faire_meta, new_rows)
        }, error = function(e) {
            faire_meta
        })

        faire_meta   <- faire_meta[order(faire_meta\$samp_name), ]

        # This assumes seqkit stats and sampsheet are ordered the same
        seqkit_stats\$sample <- sampsheet\$samp_name

        numrows         <- nrow(faire_meta)
        samp_name       <- faire_meta\$samp_name
        assay           <- rep(assay, numrows)
        pcr_plate_id    <- rep("", numrows)
        lib_id          <- fileinfo\$lib_id
        seq_run_id      <- rep("", numrows)
        lib_conc        <- rep("", numrows)
        lib_conc_unit   <- rep("", numrows)
        lib_conc_method <- rep("", numrows)
        phix_perc       <- rep("", numrows)

        if ("fw_index" %in% sampsheet) {
            mid_forward <- sampsheet\$fw_index
        } else {
            mid_forward <- rep("", numrows)
        }

        if ("rv_index" %in% sampsheet) {
            mid_reverse <- sampsheet\$rv_index
        } else {
            mid_reverse <- rep("", numrows)
        }

        filename            <- fileinfo\$filename
        filename2           <- fileinfo\$filename2
        checksum_filename   <- fileinfo\$checksum_filename
        checksum_filename2  <- fileinfo\$checksum_filename2
        associatedSequences <- rep("", numrows)
        input_read_count    <- seqkit_stats\$num_seqs

        output_readcount     <- c()
        output_otu_num       <- c()
        otu_num_tax_assigned <- c()
        for (sam in samp_name) {
            if (sam %in% colnames(phy_otu)) {
                output_readcount <- c(output_readcount, sum(phy_otu[, sam]))
                output_otu_num   <- c(output_otu_num, sum(phy_otu[, sam] > 0))
                otu_num_tax_assigned   <- c(otu_num_tax_assigned, sum(phy_tax[phy_otu[, sam] > 0, ]\$LCA != "NA"))
            } else {
                output_readcount <- c(output_readcount, 0)
                output_otu_num   <- c(output_otu_num, 0)
                otu_num_tax_assigned   <- c(otu_num_tax_assigned, 0)
            }
        }

        experimentRunMetadata <- data.frame(
            samp_name,
            assay,
            pcr_plate_id,
            lib_id,
            seq_run_id,
            lib_conc,
            lib_conc_unit,
            lib_conc_method,
            phix_perc,
            mid_forward,
            mid_reverse,
            filename,
            filename2,
            checksum_filename,
            checksum_filename2,
            associatedSequences,
            input_read_count,
            output_readcount,
            output_otu_num,
            otu_num_tax_assigned
        )

        # Write to "experimentRunMetadata" sheet
        result = tryCatch({
            writeData(wb, sheet = "experimentRunMetadata", x = experimentRunMetadata, startRow = getLastRow(wb, "experimentRunMetadata") + 1, colNames = FALSE)
        }, error = function(e) {
            addWorksheet(wb, sheetName = "experimentRunMetadata")
            writeData(wb, sheet = "experimentRunMetadata", x = experimentRunMetadata, rowNames = TRUE)
        })

        # Save workbook
        saveWorkbook(wb, paste0(${prefix}, "_faire_metadata.xlsx"), overwrite = TRUE)
    }


    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    openxlsx: ", packageVersion("openxlsx"))), "versions.yml")
    """
}
