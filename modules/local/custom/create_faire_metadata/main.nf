process CREATE_FAIRE_METADATA {
    tag "$prefix"
    label 'process_low'
    //container 'docker.io/pedrofeijao/pandas-openpyxl:v1.0'
    container 'docker.io/pawelqs/tidyverse_jsonlite_openxlsx:v1'

    input:
    tuple val(prefix), path(taxa_raw)
    path(metadata)

    output:
    path "*xlsx"       , emit: xlsx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "\"${prefix}\""
    def taxa_raw = "\"${taxa_raw}\""
    //def taxa_final = "\"${taxa_final}\""
    def metadata = "\"${metadata}\""
    """
    #!/usr/bin/env Rscript

    # Load libraries
    library(openxlsx)

    getLastRow <- function(wb, sheet) {
        nrows <- nrow(readWorkbook(wb, sheet = sheet))
        if (is.na(nrows)) return(1)
        return(nrows + 1)
    }

    # Get R version
    r_version <- R.version.string

    # Read data
    taxa_raw <- read.table(${taxa_raw}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE)
    taxa_final <- read.table(${taxa_raw}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

    # Load existing workbook
    wb <- loadWorkbook(${metadata})

    # Write to "taxaRaw" sheet
    writeData(wb, sheet = "taxaRaw", x = taxa_raw, startRow = getLastRow(wb, "taxaRaw") + 1, colNames = FALSE)

    # Write to "taxaFinal" sheet
    writeData(wb, sheet = "taxaFinal", x = taxa_final, startRow = getLastRow(wb, "taxaFinal") + 1, colNames = FALSE)

    # Save workbook
    saveWorkbook(wb, paste0(${prefix}, "_faire_metadata.xlsx"), overwrite = TRUE)

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    openxlsx: ", packageVersion("openxlsx"))), "versions.yml")
    """
}
