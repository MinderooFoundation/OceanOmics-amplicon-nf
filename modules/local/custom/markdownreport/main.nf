process MARKDOWN_REPORT {
    tag ""
    label 'process_low'

    container 'adbennett/amplicon_report:v1'

    input:
    path final_stats
    path raw_stats
    path taxa
    path pngs
    path missing

    output:
    path "*.html"      , emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def final_stats = "\"${final_stats}\""
    def raw_stats = "\"${raw_stats}\""
    def taxa = "\"${taxa}\""
    def pngs = "\"${pngs}\""
    def projectDir = "\"${projectDir}\""
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(rmarkdown))
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(png))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(scales))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(tibble))

    final_stats <- sort(c(strsplit($final_stats, " ")[[1]]))
    raw_stats   <- sort(c(strsplit($raw_stats, " ")[[1]]))
    taxa        <- sort(c(strsplit($taxa, " ")[[1]]))
    pngs        <- sort(c(strsplit($pngs, " ")[[1]]))

    if (length(raw_stats) > 0) {
        concat_raw_df     <- data.frame()

        for (i in raw_stats) {
            curr_table    <- read_table(i)
            concat_raw_df <- rbind(concat_raw_df, curr_table)
        }

        concat_raw_df     <- concat_raw_df[!duplicated(concat_raw_df), ]
        write.table(concat_raw_df, file = "concat_raw_stats.tsv", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

    if (length(final_stats) > 0) {
        concat_sample_df     <- data.frame()

        for (i in final_stats) {
            curr_table       <- read_table(i)
            concat_sample_df <- rbind(concat_sample_df, curr_table)
        }

        write.table(concat_sample_df, file = "concat_sample_stats.tsv", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }

    # Copy the R markdown scripts to this work directory to avoid issues with R markdown changing the work directory
    file.copy(paste0(${projectDir}, "/bin/amplicon_report.Rmd"), ".")

    # Render the html file
    render("amplicon_report.Rmd", output_file = "amplicon_report_mqc.html")

    # Version information
    writeLines(c("\\"${task.process}\\":",
                paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                paste0("    rmarkdown: ", packageVersion("rmarkdown")),
                paste0("    tidyverse: ", packageVersion("tidyverse")),
                paste0("    png: ", packageVersion("png")),
                paste0("    ggplot2: ", packageVersion("ggplot2")),
                paste0("    scales: ", packageVersion("scales")),
                paste0("    reshape2: ", packageVersion("reshape2")),
                paste0("    tibble: ", packageVersion("tibble"))),
                "versions.yml")
    """
}
