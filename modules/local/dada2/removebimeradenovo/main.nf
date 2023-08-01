process DADA2_REMOVEBIMERADENOVO {
    tag "$ids"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids
    path seq_table

    output:
    path("seq_table_nochim.RData"), emit: seq_table_nochim
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seq_table = "\"${seq_table}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    seq_table <- readRDS($seq_table)

    # Remove Chimeras
    seq_table_nochim <- removeBimeraDenovo(seq_table, 
                                           method = "pooled", 
                                           multithread = ${task.cpus},
                                           verbose = TRUE)

    if (length(colnames(seq_table_nochim)) < 2) {
        print(paste0("Only ", length(colnames(seq_table_nochim)), " ASVs left after removing chimera, at least 2 ASVs are needed; a larger dataset may be needed. You can also try changing DADA2 options ('--pooled', '--min_overlap', or '--max_mismatch')"))
        stop()
    }
  
    # Save the result
    save(seq_table_nochim, file = "seq_table_nochim.RData")

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}