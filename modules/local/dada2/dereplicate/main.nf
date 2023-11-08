process DADA2_DEREPLICATE {
    tag "$ids"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids
    val single_end
    path fq_files

    output:
    path("dereplicated.RData"), emit: dereplicate
    path "versions.yml"       , emit: versions

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

    # Dereplication
    if (! single_end) {
        fw_list <- grep("*R1_trimmed.fq.gz", fq_list, value = TRUE)
        rv_list <- grep("*R2_trimmed.fq.gz", fq_list, value = TRUE)
        derep_forward        <- derepFastq(fw_list, verbose = TRUE)
        derep_reverse        <- derepFastq(rv_list, verbose = TRUE)
        names(derep_forward) <- id_list
        names(derep_reverse) <- id_list

    } else {
        fw_list <- grep("*trimmed.fq.gz", fq_list, value = TRUE)
        derep_forward        <- derepFastq(fw_list, verbose = TRUE)
        names(derep_forward) <- id_list
    }

    # Save dereplication
    if (single_end) {
        save(derep_forward, file = paste0("dereplicated.RData"))
    } else {
        save(derep_forward, derep_reverse, file = paste0("dereplicated.RData"))
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}
