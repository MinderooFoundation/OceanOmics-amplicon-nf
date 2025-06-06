process CREATE_DEMUX_DEPENDENCIES {
    tag "$index_file"
    label 'process_medium'
    container 'rocker/tidyverse:4.3.0'

    input:
    path index_file

    output:
    path "fw.fa"             , emit: fw_index
    path "rv.fa"             , emit: rv_index
    path "rename_pattern.txt", emit: sample_rename
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def index_file = "\"$index_file\""
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(dplyr))

    barcodes <- read_csv(${index_file}, col_types=cols(.default='c'))

    if ("discarded" %in% colnames(barcodes)) {
        barcodes <- barcodes[barcodes\$discarded == FALSE, ]
    }

    # If there is no rv_index then it is single end
    if (any(is.na(barcodes\$rv_index))) {
        barcodes <- barcodes %>%
            group_by(fw_index) %>%
            mutate(fw_no = paste0("F", cur_group_id())) %>%
            ungroup()

        barcodes_fw <- barcodes %>%
            distinct(fw_index, .keep_all = TRUE, .keep_last = FALSE)

        names   <- barcodes_fw\$fw_no
        indexes <- barcodes_fw\$fw_index

        cat(
            paste(
                paste0(">", names),
                indexes,
                sep="\n"
            ),
            sep = "\n",
            file = "fw.fa"
        )
        cat(
            paste(
                paste0(">", names),
                indexes,
                sep="\n"
            ),
            sep = "\n",
            file = "rv.fa"
        )

        cat(paste0(barcodes\$fw_no, ".R1.fq.gz ", barcodes\$samp_name , ".1.fq.gz"),
            sep="\n",
            file="rename_pattern.txt")

    } else {
        # paired end
        barcodes <- barcodes %>%
            group_by(fw_index) %>%
            mutate(fw_no = paste0("F", cur_group_id())) %>%
            ungroup()
        barcodes <- barcodes %>%
            group_by(rv_index) %>%
            mutate(rv_no = paste0("R", cur_group_id())) %>%
            ungroup()

        barcodes_fw <- barcodes %>%
            distinct(fw_index, .keep_all = TRUE, .keep_last = FALSE)
        barcodes_rv <- barcodes %>%
            distinct(rv_index, .keep_all = TRUE, .keep_last = FALSE)

        fw_names   <- barcodes_fw\$fw_no
        rv_names   <- barcodes_rv\$rv_no
        names      <- c(fw_names, rv_names)
        fw_indexes <- barcodes_fw\$fw_index
        rv_indexes <- barcodes_rv\$rv_index
        indexes    <- c(fw_indexes, rv_indexes)

        # This will generate a .fa file that searches for both the Fw and the Rv file in R1; whilst keeping the same sample name.
        cat(
            paste(
                paste0(">", names),
                indexes,
                sep="\n"
            ),
            sep = "\n",
            file = "fw.fa"
        )

        cat(
            paste(
                paste0(">",names),
                indexes,
                sep="\n"
            ),
            sep = "\n",
            file = "rv.fa"
        )

        # This section creates a file to rename the demultiplexed files to reflect the sample name, including the assay
        cat(paste0(barcodes\$fw_no, "-", barcodes\$rv_no, ".R[12].fq.gz ", barcodes\$samp_name , "_forward.#1.fq.gz"),
            paste0(barcodes\$rv_no, "-", barcodes\$fw_no, ".R[12].fq.gz ", barcodes\$samp_name , "_reverse.#1.fq.gz"),
            sep = "\n",
            file = "rename_pattern.txt")
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    tidyverse: ", packageVersion("tidyverse")),paste0("    dplyr: ", packageVersion("dplyr"))), "versions.yml")
    """
}
