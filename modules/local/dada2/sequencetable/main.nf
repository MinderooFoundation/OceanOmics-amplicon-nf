process DADA2_SEQUENCETABLE {
    tag "$ids"
    label 'process_medium'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids 
    val single_end
    path sample_inference
    path merged

    output:
    path "seq_tab.rds" , emit: seq_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def single_end = "\"${single_end}\""
    def sample_inference = "\"${sample_inference}\""
    def merged = "\"${merged}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    single_end             <- as.logical($single_end)

    if (single_end) {
        load($sample_inference)
        
        seq_table <- makeSequenceTable(dada_forward)

        saveRDS(seq_table, file = "seq_tab.rds")
        
    } else {
        load($merged)

        seq_table <- makeSequenceTable(mergers)

        saveRDS(seq_table, file = "seq_tab.rds")
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}