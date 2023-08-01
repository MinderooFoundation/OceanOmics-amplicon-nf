process DADA2_SAMPLEINFERENCE {
    tag "$ids"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids 
    val single_end
    path dereplicate
    path error_rates

    output:
    path("sample_inference.RData"), emit: sample_inference
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def single_end = "\"${single_end}\""
    def dereplicate = "\"${dereplicate}\""
    def error_rates = "\"${error_rates}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    load($dereplicate)
    load($error_rates)
    single_end             <- as.logical($single_end)
    
    # Sample inference
    dada_forward <- dada(derep_forward, 
                         err = errors_forward, 
                         $args, 
                         multithread = ${task.cpus},
                         verbose = TRUE)

    if (! single_end) {
    dada_reverse <- dada(derep_reverse, 
                         err = errors_reverse, 
                         $args, 
                         multithread = ${task.cpus},
                         verbose = TRUE)
    }
    
    # Save sample inferences
    if (single_end) {
        save(dada_forward, file = "sample_inference.RData")
    } else {
        save(dada_forward, dada_reverse, file = "sample_inference.RData")
    }
  
    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}