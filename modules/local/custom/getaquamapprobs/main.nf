process GET_AQUAMAP_PROBS {
    tag "$samplesheet"
    label 'process_medium'
    container 'adbennett/phyloseq_ncdf4_raster:v1.0'

    input:
    tuple val(prefix), path(phyloseq)
    tuple val(prefix), path(nc_files)

    output:
    tuple val(prefix), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    library(ncdf4)
    library(phyloseq)

    if (file.exists("dummy.nc")) {
        file.create(paste0("no_aquamaps_", "$phyloseq", ".csv"))
    } else {
        phyloseq         <- readRDS("$phyloseq")
        species          <- unique(phyloseq@tax_table@.Data[, "species"])
        species          <- species[!is.na(species)]
        species          <- species[species != "dropped"]
        samples          <- rownames(phyloseq@sam_data)
        out_df           <- data.frame(matrix(ncol = length(samples), nrow = length(species)))
        colnames(out_df) <- samples
        rownames(out_df) <- species

        if ("latitude" %in% colnames(phyloseq@sam_data) & "longitude" %in% colnames(phyloseq@sam_data)) {

            for (spec in species) {
                spec <- gsub(" ", "_", spec)

                if (file.exists(paste0(spec, ".nc"))) {
                    nc              <- nc_open(paste0(spec, ".nc"))
                    probs           <- data.frame(ncvar_get(nc, varid = "probability"))
                    lats            <- ncvar_get(nc, varid = "latitude")
                    longs           <- ncvar_get(nc, varid = "longitude")
                    colnames(probs) <- lats
                    rownames(probs) <- longs

                    for (sam in samples) {
                        sample_lat  <- data.frame(phyloseq@sam_data)[sam, "latitude"]
                        sample_long <- data.frame(phyloseq@sam_data)[sam, "longitude"]
                        if (! is.na(sample_lat) & ! is.na(sample_long)) {
                            prob_na     <- FALSE

                            # Round lat/long to the closest edge if it's already close to the edge
                            # Assign NA to prob if the lat/long is too far outside the range of lats/longs
                            if (sample_lat < min(lats) & sample_lat >= min(lats) - 0.5) {
                                sample_lat <- min(lats)

                            } else if (sample_lat > max(lats) & sample_lat <= max(lats) + 0.5) {
                                sample_lat <- max(lats)

                            } else if (sample_lat < min(lats) | sample_lat > max(lats)) {
                                prob    <- NA
                                prob_na <- TRUE
                            }

                            if (sample_long < min(longs) & sample_long >= min(longs) - 0.5 & prob_na != TRUE) {
                                sample_long <- min(longs)

                            } else if (sample_long > max(longs) & sample_long <= max(longs) + 0.5) {
                                sample_long <- max(longs)

                            } else if (sample_long < min(longs) | sample_long > max(longs)) {
                                prob    <- NA
                                prob_na <- TRUE
                            }

                            # Round lat/long to the closest lat/long
                            if (prob_na != TRUE) {
                                if (! sample_lat %in% lats) {
                                    sample_lat <- lats[which.min(abs(lats - sample_lat))]
                                }

                                if (! sample_long %in% longs) {
                                    sample_long <- longs[which.min(abs(longs - sample_long))]
                                }

                                prob <- probs[toString(sample_long), toString(sample_lat)]
                            }

                            spec <- gsub("_", " ", spec)
                            out_df[spec, sam] <- prob
                        } else {
                            spec <- gsub("_", " ", spec)
                            out_df[spec, sam] <- NA
                        }
                    }
                } else {
                    out_df[spec, ] <- "Species not in aquamaps"
                }
            }

            # Output csv
            write.csv(out_df, paste0("aquamaps_", "$phyloseq", ".csv"))

        } else {
            file.create(paste0("no_aquamaps_", "$phyloseq", ".csv"))
        }
    }
    """
}
