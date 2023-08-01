process DADA2_PLOTQUALITYPROFILE {
    tag "$meta"
    label 'process_single'
    container 'quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0'

    input:
    tuple val(meta), path(reads)
    val stage

    output:
    tuple val(meta), path("*.png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "\"${meta.id}\""
    def fq_files  = meta.single_end ? "\"${reads}\"" : "\"${reads[0]}, ${reads[1]}\""
    def single_end = "\"${meta.single_end}\""
    def stage = "\"${stage}\""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    single_end             <- as.logical($single_end)
    
    if (single_end) {
        # Create plot
        quality_profile    <- plotQualityProfile($fq_files)
        
        # Save plot
        png(paste0($prefix, "_", $stage, ".png"))
        print(quality_profile)
        dev.off()

    } else {
        # Create plots
        fq_list            <- as.list(strsplit($fq_files, ", ")[[1]])
        quality_profile_fw <- plotQualityProfile(fq_list[1])
        quality_profile_rv <- plotQualityProfile(fq_list[2])

        # Save plots
        png(paste0($prefix, "_", $stage, "_R1.png"))
        print(quality_profile_fw)
        dev.off()
        png(paste0($prefix, "_", $stage, "_R2.png"))
        print(quality_profile_rv)
        dev.off()
    }

    # Version information
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2"))), "versions.yml")
    """
}