process DADA2_PLOTERRORS {
    tag "$ids"
    label 'process_single'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    val ids
    val single_end
    path error_rates

    output:
    path "*.png"       , emit: png
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def single_end = "\"${single_end}\""
    def error_rates = "\"${error_rates}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    load($error_rates)
    single_end             <- as.logical($single_end)

    # Create plot
    errors_plot_fw <- plotErrors(errors_forward, nominalQ = TRUE)

    # Save plot
    png("errors_plot_fw.png")
    print(errors_plot_fw)
    dev.off()

    if (! single_end) {
        # Create plot
        errors_plot_rv <- plotErrors(errors_reverse, nominalQ = TRUE)

        # Save plot
        png("errors_plot_rv.png")
        print(errors_plot_rv)
        dev.off()
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}
