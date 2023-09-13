process ZOTUTRACKREADS {
    tag ""
    label 'process_medium'

    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val prefix
    path merged_counts
    path derep_counts
    path denoise_counts

    output:
    path "*.txt"        , emit: txt
    path "*.png"        , emit: png
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "\"${prefix}\""
    def merged_counts = "\"${merged_counts}\""
    def derep_counts = "\"${derep_counts}\""
    def denoise_counts = "\"${denoise_counts}\""
    """
    #!/usr/bin/env Rscript

    # Read the merged count file
    merged_file  <- readLines(${merged_counts})
    sample_names <- unique(merged_file[seq(1, length(merged_file), by = 2)])
    merged_nums  <- as.numeric(merged_file[seq(2, length(merged_file), by = 2)])
    merged_nums  <- append(merged_nums, NA)

    # Read the dereplicated count file
    derep_file   <- readLines(${derep_counts})
    d_samp_names <- derep_file[seq(1, length(derep_file), by = 2)]
    derep_nums   <- as.numeric(derep_file[seq(2, length(derep_file), by = 2)])
    derep_nums   <- append(derep_nums, NA)

    missing_samples <- setdiff(sample_names, d_samp_names)
    for (missing_sample in missing_samples) {
        derep_nums <- append(derep_nums, 0)
    }

    # Read the denoised count file
    denoise_file <- readLines(${denoise_counts})
    denoise_nums <- as.numeric(denoise_file[seq(2, length(denoise_file), by = 2)])

    # Add 'zotus' to the sample names for the data frame
    all_sample_names <- append(sample_names, "zotus")

    # Create an empty data frame of the correct size
    track_df <- data.frame(merged = rep(NA, length(all_sample_names)),
                            dereplicated = rep(NA, length(all_sample_names)),
                            denoised = rep(NA, length(all_sample_names)))
    rownames(track_df) <- all_sample_names

    # Add the count numbers to the data frame
    track_df\$merged              <- merged_nums
    track_df\$dereplicated        <- derep_nums
    track_df["zotus", "denoised"] <- sum(denoise_nums)

    #track reads
    write.table(track_df, file = paste0(${prefix}, "_track_reads.txt"))

    # Modify data frame for box plot
    samps                <- row.names(track_df)
    track_df\$samps      <- samps
    track_df_long        <- reshape(track_df, varying = list(c("merged", "dereplicated", "denoised")), direction = "long", timevar = "stage", times = c("merged", "dereplicated", "denoised"), v.names = "reads")
    track_df_long\$stage <- factor(track_df_long\$stage, levels = unique(track_df_long\$stage))

    # Save plot
    png(paste0(${prefix}, "_track_reads.png"), width=800)
    boxplot(reads ~ stage, data = track_df_long, col = "lightblue", boxwex = 0.5, xlab="Sequencing stage", ylab="Number of reads")
    dev.off()

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = "."))), "versions.yml")
    """
}
