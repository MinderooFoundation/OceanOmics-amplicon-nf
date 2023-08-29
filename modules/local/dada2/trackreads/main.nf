process DADA2_TRACKREADS {
    tag ""
    label 'process_medium'

    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val prefix
    val single_end
    path sample_inference
    path seq_table_nochim
    path merged
    path ft_out

    output:
    path "*.txt"       , emit: txt
    path "*.png"       , emit: png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "\"${prefix}\""
    def single_end = "\"${single_end}\""
    def sample_inference = "\"${sample_inference}\""
    def seq_table_nochim = "\"${seq_table_nochim}\""
    def merged = "\"${merged}\""
    def ft_out = "\"${ft_out}\""
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(dada2))

    load($sample_inference)
    load($seq_table_nochim)

    # Collect all filter and trim out csv files into one data frame
    csv_files     <- sort(c(strsplit($ft_out, " ")[[1]]))
    out           <- do.call(rbind, lapply(csv_files, read.csv))
    rownames(out) <- out[,1]
    out[,1]       <- NULL

    single_end    <- as.logical($single_end)

    # Helper function to get sum of uniques
    getN <- function(x) sum(getUniques(x))

    if (single_end) {
        #track reads
        track_df <- cbind(out, sapply(dada_forward, getN), rowSums(seq_table_nochim))
        track_df <- as.data.frame(track_df)
        colnames(track_df) <- c("input", "filtered", "denoised", "nonchim")
        rownames(track_df) <- rownames(out)

        write.table(track_df, file = paste0(${prefix}, "_track_reads.txt"))

        # Modify data frame for box plot
        samps                <- row.names(track_df)
        track_df\$samps      <- samps
        track_df_long        <- reshape(track_df, varying = list(c("input", "filtered", "denoised", "nonchim")), direction = "long", timevar = "stage", times = c("input", "filtered", "denoised", "nonchim"), v.names = "reads")
        track_df_long\$stage <- factor(track_df_long\$stage, levels = unique(track_df_long\$stage))

        # Save plot
        png(paste0(${prefix}, "_track_reads.png"), width=800)
        boxplot(reads ~ stage, data = track_df_long, col = "lightblue", boxwex = 0.5, xlab="Sequencing stage", ylab="Number of reads")
        dev.off()

    # Paired end
    } else {
        load($merged)

        # Track reads
        track_df <- cbind(out, sapply(dada_forward, getN), sapply(dada_reverse, getN), sapply(mergers, getN), rowSums(seq_table_nochim))
        track_df <- as.data.frame(track_df)
        colnames(track_df) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        rownames(track_df) <- rownames(out)

        write.table(track_df, file = paste0(${prefix}, "_track_reads.txt"))

        # Modify data frame for box plot
        samps                <- row.names(track_df)
        track_df\$samps      <- samps
        track_df_long        <- reshape(track_df, varying = list(c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")), direction = "long", timevar = "stage", times = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"), v.names = "reads")
        track_df_long\$stage <- factor(track_df_long\$stage, levels = unique(track_df_long\$stage))

        # Save plot
        png(paste0(${prefix}, "_track_reads.png"), width=800)
        boxplot(reads ~ stage, data = track_df_long, col = "lightblue", boxwex = 0.5, xlab="Sequencing stage", ylab="Number of reads")
        dev.off()
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = "."))), "versions.yml")
    """
}
