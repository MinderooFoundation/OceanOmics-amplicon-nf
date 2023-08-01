process DADA2_FINALOUTPUTS {
    tag "$ids"
    label 'process_low'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val prefix
    val ids
    path seq_table_nochim

    output:
    tuple val(prefix), path("asv_table.csv")      , emit: asv_csv
    tuple val(prefix), path("asv.fa")             , emit: asv_fasta
    tuple val(prefix), path("asv_final_table.tsv"), emit: asv_tsv
    tuple val(prefix), path("lca_input.tsv")      , emit: lca_input_tsv
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def seq_table_nochim = "\"${seq_table_nochim}\""
    """
    #!/usr/bin/env Rscript

    #suppressPackageStartupMessages(library(dada2))

    load($seq_table_nochim)

    if (is.null(rownames(seq_table_nochim))) {
        cat("Error: Couldn't create DADA2 final outputs. At least two samples are needed.")
        q(status = 1)
    }

    asv_seqs    <- colnames(seq_table_nochim)
    asv_headers <- vector(dim(seq_table_nochim)[2], mode="character")
  
    for (i in 1:dim(seq_table_nochim)[2]) {
        asv_headers[i] <- paste(">ASV", i, sep="_")
    }
  
    # Format the asv table
    asv_table            <- seq_table_nochim
    colnames(asv_table)  <- asv_headers
    IDs                  <- rownames(asv_table)
    asv_table            <- as.data.frame(asv_table)
    colnames(asv_table)  <- asv_headers
    asv_table\$sample_id <- IDs
    asv_table            <- asv_table[, c("sample_id", asv_headers)]
  
    # Save the asv table
    write.table(asv_table, file = "asv_table.csv", sep = ",", row.names = FALSE, quote = FALSE)

    # Create and save asv fasta file
    asv_fasta <- c(rbind(asv_headers, asv_seqs))
    write(asv_fasta, file = "asv.fa")

    # Now we need to format the table for the final outputs
    colnames(seq_table_nochim) <- asv_headers
  
    # We need to transpose, so the rows are the sequences, whereas the columns are the IDs
    asv_for_lca <- as.data.frame(t(seq_table_nochim))
  
    # Add asv column and move column to left of table
    asv_for_lca\$ASV <- rownames(asv_for_lca)
    asv_for_lca <- asv_for_lca[, c(ncol(asv_for_lca), 1:(ncol(asv_for_lca)-1))]

    # Remove '>' and add sequences
    asv_for_lca[, 1] <- gsub(">", "", asv_for_lca[, 1])
    asv_for_lca\$ASV_sequence <- asv_seqs

    # Save the final asv table
    write.table(asv_for_lca, file = "asv_final_table.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

    # Format table for lca input
    lca_input <- asv_for_lca
    names(lca_input)[1] <- "#ID"
    lca_input <- lca_input[, -ncol(asv_for_lca)]
    
    # Save the lca input table
    write.table(lca_input, file = "lca_input.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}