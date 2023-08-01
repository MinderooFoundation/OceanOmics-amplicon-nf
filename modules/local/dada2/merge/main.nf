process DADA2_MERGE {
    tag "$ids"
    label 'process_medium'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids 
    val single_end
    path sample_inference
    path fq_files

    output:
    path("merged*")    , emit: merged
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def single_end = "\"${single_end}\""
    def sample_inference = "\"${sample_inference}\""
    def fq_files = "\"${fq_files}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    single_end             <- as.logical($single_end)

    if (single_end) {
        fileConn<-file("merged.txt")
        writeLines(c("The data is single end", "fq files were not merged"), fileConn)
        close(fileConn)
        
    } else {
        load($sample_inference)

        fq_list <- sort(c(strsplit($fq_files, " ")[[1]]))
        fw_list <- grep("*R1_trimmed.fq.gz", fq_list, value = TRUE)
        rv_list <- grep("*R2_trimmed.fq.gz", fq_list, value = TRUE)

        mergers <- mergePairs(dada_forward, 
                              fw_list, 
                              dada_reverse, 
                              rv_list,
                              $args, 
                              verbose=TRUE)

        save(mergers, file = "merged.RData")
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}