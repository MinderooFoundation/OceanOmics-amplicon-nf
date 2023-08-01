process DADA2_SEQUENCEDISTRIBUTION {
    tag "$ids"
    label 'process_single'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids
    path seq_table

    output:
    path("*.png")      , emit: png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seq_table = "\"${seq_table}\""
    """
    #!/usr/bin/env Rscript

    seq_table <- readRDS($seq_table)

    seq_lengths <- nchar(colnames(seq_table))
  
    # Create and save plot
    png("ASV_seq_distribution.png")
    hist(seq_lengths, breaks = 100, xlab = "Sequence length (bp)", ylab = "Number of reads", main = "", col = "lightblue")
    dev.off()

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}