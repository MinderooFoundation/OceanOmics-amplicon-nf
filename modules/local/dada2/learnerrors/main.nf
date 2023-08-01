process DADA2_LEARNERRORS {
    tag "$ids"
    label 'process_medium'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids
    val single_end
    path fq_files 

    output:
    path("error_rates.RData"), emit: error_rates
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = "\"${ids}\"" 
    def single_end = "\"${single_end}\""
    def fq_files = "\"${fq_files}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    single_end             <- as.logical($single_end)

    # remove '[' and ']' from input string, then split into lists
    id_list <- gsub('^.|.\$', '', $ids)
    id_list <- sort(c(strsplit(id_list, ", ")[[1]]))
    fq_list <- sort(c(strsplit($fq_files, " ")[[1]]))

    # Learn error rates
    if (! single_end) {
        fw_list <- grep("*R1_trimmed.fq.gz", fq_list, value = TRUE)
        rv_list <- grep("*R2_trimmed.fq.gz", fq_list, value = TRUE)
        errors_forward <- learnErrors(fw_list, multithread = $task.cpus)
        errors_reverse <- learnErrors(rv_list, multithread = $task.cpus)
    } else {
        fw_list <- grep("*trimmed.fq.gz", fq_list, value = TRUE)
        errors_forward <- learnErrors(fw_list, multithread = $task.cpus)
    }
    
    # Save error rates
    if (single_end) {
        save(errors_forward, file = paste0("error_rates.RData"))
    } else {
        save(errors_forward, errors_reverse, file = paste0("error_rates.RData"))
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}