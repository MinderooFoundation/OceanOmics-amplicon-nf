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

    # If there is no rv_index then it is single end
    if (any(is.na(barcodes\$rv_index))) {
        barcodes <- barcodes %>%
            group_by(fw_index) %>%
            mutate(fw_no = paste0("F", cur_group_id())) %>%
            ungroup()

        cat(paste(paste0(">", barcodes\$fw_no),
            barcodes\$fw_index, sep="\n"),
            sep = "\n", 
            file = "fw.fa")
        cat(paste(paste0(">", barcodes\$fw_no),
            barcodes\$fw_index, sep="\n"),
            sep = "\n", 
            file = "rv.fa")

        cat(paste0(barcodes\$fw_no, ".R1.fq.gz ", barcodes\$sample , ".1.fq.gz"),
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
    
        # This will generate a .fa file that searches for both the Fw and the Rv file in R1; whilst keeping the same sample name.
        cat(paste(paste0(">", barcodes\$fw_no),
            barcodes\$fw_index,
            paste0(">", barcodes\$rv_no),
            barcodes\$rv_index, sep="\n"), 
            sep = "\n", 
            file = "fw.fa")

        cat(paste(paste0(">",barcodes\$rv_no),
            barcodes\$rv_index,
            paste0(">",barcodes\$fw_no),
            barcodes\$fw_index, sep="\n"), 
            sep = "\n", 
            file = "rv.fa")

        # This section creates a file to rename the demultiplexed files to reflect the sample name, including the assay
        cat(paste0(barcodes\$fw_no, "-", barcodes\$rv_no, ".R[12].fq.gz ", barcodes\$sample , "_forward.#1.fq.gz"),
            paste0(barcodes\$rv_no, "-", barcodes\$fw_no, ".R[12].fq.gz ", barcodes\$sample , "_reverse.#1.fq.gz"),
            sep = "\n",
            file = "rename_pattern.txt")
    }
 
    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    tidyverse: ", packageVersion("tidyverse")),paste0("    dplyr: ", packageVersion("dplyr"))), "versions.yml")
    """
}