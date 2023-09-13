process DADA2_FILTERANDTRIM {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*trimmed.fq.gz"), emit: reads
    path "*.csv"                           , emit: csv
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "\"${meta.id}\""
    def fq_files  = meta.single_end ? "\"${reads}\"" : "\"${reads[0]}, ${reads[1]}\""
    def single_end = "\"${meta.single_end}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    single_end             <- as.logical($single_end)

    # Single end
    if (single_end) {
        ft_out <- filterAndTrim($fq_files,
                                paste0($prefix, "_trimmed.fq.gz"),
                                multithread=$task.cpus)

    # Paired end
    } else {
        fq_list            <- c((strsplit($fq_files, ", ")[[1]]))
        out_list           <- c(paste0($prefix, "_R1_trimmed.fq.gz"), paste0($prefix, "_R2_trimmed.fq.gz"))

        ft_out <- filterAndTrim(fq_list[1],
                                out_list[1],
                                fq_list[2],
                                out_list[2],
                                multithread=$task.cpus)
    }

    writeLines(c("sample,reads.in,reads.out", paste0($prefix, ",", ft_out[1,1], ",", ft_out[1,2])), paste0($prefix, "_filter_trim_out.csv"))

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}
