process CREATE_FAIRE_METADATA {
    tag "$prefix"
    label 'process_low'
    //container 'docker.io/pedrofeijao/pandas-openpyxl:v1.0'
    container 'docker.io/pawelqs/tidyverse_jsonlite_openxlsx:v1'

    input:
    tuple val(prefix), path(taxa_raw), path(taxa_final), path(otu_raw), path(otu_final)
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
    def taxa_final = "\"${taxa_final}\""
    def otu_raw = "\"${otu_raw}\""
    def otu_final = "\"${otu_final}\""
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
    taxa_raw <- read.table(${taxa_raw}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")
    taxa_final <- read.table(${taxa_final}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")
    otu_raw <- read.table(${otu_raw}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")
    otu_final <- read.table(${otu_final}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")

    # Load existing workbook
    wb <- loadWorkbook(${metadata})

    # Write to "taxaRaw" sheet
    result = tryCatch({
        writeData(wb, sheet = "taxaRaw", x = taxa_raw, startRow = getLastRow(wb, "taxaRaw") + 1, colNames = FALSE)
    }, error = function(e) {
        addWorksheet(wb, sheetName = "taxaRaw")
        writeData(wb, sheet = "taxaRaw", x = taxa_raw, rowNames = TRUE)
    })

    # Write to "taxaFinal" sheet
    result = tryCatch({
        writeData(wb, sheet = "taxaFinal", x = taxa_final, startRow = getLastRow(wb, "taxaFinal") + 1, colNames = FALSE)
    }, error = function(e) {
        addWorksheet(wb, sheetName = "taxaFinal")
        writeData(wb, sheet = "taxaFinal", x = taxa_final, rowNames = TRUE)
    })

    # Write to "otuRaw" sheet
    result = tryCatch({
        writeData(wb, sheet = "otuRaw", x = otu_raw, startRow = getLastRow(wb, "otuRaw") + 1, colNames = FALSE)
    }, error = function(e) {
        addWorksheet(wb, sheetName = "otuRaw")
        writeData(wb, sheet = "otuRaw", x = otu_raw, rowNames = TRUE)
    })

    # Write to "otuFinal" sheet
    result = tryCatch({
        writeData(wb, sheet = "otuFinal", x = otu_final, startRow = getLastRow(wb, "otuFinal") + 1, colNames = FALSE)
    }, error = function(e) {
        addWorksheet(wb, sheetName = "otuFinal")
        writeData(wb, sheet = "otuFinal", x = otu_final, rowNames = TRUE)
    })

    # Save workbook
    saveWorkbook(wb, paste0(${prefix}, "_final_faire_metadata.xlsx"), overwrite = TRUE)

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    openxlsx: ", packageVersion("openxlsx"))), "versions.yml")
    """
}
