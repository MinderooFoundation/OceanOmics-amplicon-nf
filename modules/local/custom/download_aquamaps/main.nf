process DOWNLOAD_AQUAMAPS {
    tag "$samplesheet"
    label 'process_medium'
    container 'adbennett/phyloseq_and_tree:v2'

    input:
    tuple val(prefix), path(phyloseq)

    output:
    tuple val(prefix), path("*.nc"), emit: nc_files

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    phyloseq    <- readRDS("$phyloseq")

    if ("latitude" %in% colnames(phyloseq@sam_data) & "longitude" %in% colnames(phyloseq@sam_data)) {
        spec_to_get <- unique(phyloseq@tax_table@.Data[, "species"])
        spec_to_get <- spec_to_get[!is.na(spec_to_get)]
        spec_to_get <- spec_to_get[spec_to_get != "dropped"]
        urls        <- c()
        destfiles   <- c()

        for(i in spec_to_get){
            i        <- gsub(" ", "_", i)
            url      <- paste0("https://thredds.d4science.org/thredds/fileServer/public/netcdf/AquaMaps_11_2019/", i, ".nc")
            destfile <- paste0("./", i, ".nc")

            if(! file.exists(paste0("./", i, ".nc")) ) {
                urls      <- c(urls, url)
                destfiles <- c(destfiles, destfile)
            }
        }

        res       <- curl::multi_download(urls, destfiles)
        delete_us <- res\$destfile[res\$status_code == "404"]
        file.remove(delete_us)
    }

    # Create a fake file if .no nc files were downloaded
    # This is to avoid a Nextflow 'Missing output file(s) error'
    if (length(list.files(pattern = "*.nc")) == 0) {
        file.create("dummy.nc")
    }
    """
}
