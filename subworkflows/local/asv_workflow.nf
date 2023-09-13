/*
 * ASV Creation with DADA2
 */

include { DADA2_DEREPLICATE                 } from '../../modules/local/dada2/dereplicate/main'
include { DADA2_FILTERANDTRIM               } from '../../modules/local/dada2/filterandtrim/main'
include { DADA2_LEARNERRORS                 } from '../../modules/local/dada2/learnerrors/main'
include { DADA2_MERGE                       } from '../../modules/local/dada2/merge/main'
include { DADA2_PLOTERRORS                  } from '../../modules/local/dada2/ploterrors/main'
include { DADA2_PLOTQUALITYPROFILE as  \
            DADA2_PLOTQUALITYPROFILERAW;
            DADA2_PLOTQUALITYPROFILE as  \
            DADA2_PLOTQUALITYPROFILETRIMMED } from '../../modules/local/dada2/plotqualityprofile/main'
include { DADA2_REMOVEBIMERADENOVO          } from '../../modules/local/dada2/removebimeradenovo/main'
include { DADA2_SAMPLEINFERENCE             } from '../../modules/local/dada2/sampleinference/main'
include { DADA2_SEQUENCEDISTRIBUTION        } from '../../modules/local/dada2/sequencedistribution/main'
include { DADA2_SEQUENCETABLE               } from '../../modules/local/dada2/sequencetable/main'
include { DADA2_TRACKREADS                  } from '../../modules/local/dada2/trackreads/main'
include { DADA2_FINALOUTPUTS                } from '../../modules/local/dada2/finaloutputs/main'

workflow ASV_WORKFLOW {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_prefix   = "asv"

    // Plot quality of raw reads
    DADA2_PLOTQUALITYPROFILERAW (
        ch_reads,
        "raw"
    )
    ch_versions = ch_versions.mix(DADA2_PLOTQUALITYPROFILERAW.out.versions.first())

    // Filter reads and trim if necessary
    DADA2_FILTERANDTRIM (
        ch_reads
    )
    ch_versions = ch_versions.mix(DADA2_FILTERANDTRIM.out.versions.first())

    // Plot quality of filtered/trimmed reads
    DADA2_PLOTQUALITYPROFILETRIMMED (
        DADA2_FILTERANDTRIM.out.reads,
        "filtered"
    )
    ch_versions = ch_versions.mix(DADA2_PLOTQUALITYPROFILETRIMMED.out.versions.first())

    // We are collecting the samples together because the rest of the DADA2 workflow needs the samples collected
    ch_ids           = DADA2_FILTERANDTRIM.out.reads.map { it[0].id }.collect()
    ch_single_end    = DADA2_FILTERANDTRIM.out.reads.map { it[0].single_end }.distinct()
    ch_reads_filt    = DADA2_FILTERANDTRIM.out.reads.map { it[1] }.collect()
    ch_filt_trim_csv = DADA2_FILTERANDTRIM.out.csv.collect()

    // Use all samples to learn error rates
    DADA2_LEARNERRORS (
        ch_ids,
        ch_single_end,
        ch_reads_filt
    )
    ch_versions = ch_versions.mix(DADA2_LEARNERRORS.out.versions)

    // Plot the error rates
    DADA2_PLOTERRORS (
        ch_ids,
        ch_single_end,
        DADA2_LEARNERRORS.out.error_rates
    )
    ch_versions = ch_versions.mix(DADA2_PLOTERRORS.out.versions)

    // Dereplicate
    DADA2_DEREPLICATE (
        ch_ids,
        ch_single_end,
        ch_reads_filt
    )
    ch_versions = ch_versions.mix(DADA2_DEREPLICATE.out.versions)

    // Sample Inference
    DADA2_SAMPLEINFERENCE (
        ch_ids,
        ch_single_end,
        DADA2_DEREPLICATE.out.dereplicate,
        DADA2_LEARNERRORS.out.error_rates
    )
    ch_versions = ch_versions.mix(DADA2_SAMPLEINFERENCE.out.versions)

    // Merge reads
    DADA2_MERGE (
        ch_ids,
        ch_single_end,
        DADA2_SAMPLEINFERENCE.out.sample_inference,
        ch_reads_filt
    )
    ch_versions = ch_versions.mix(DADA2_MERGE.out.versions)

    // Create sequence table
    DADA2_SEQUENCETABLE (
        ch_ids,
        ch_single_end,
        DADA2_SAMPLEINFERENCE.out.sample_inference,
        DADA2_MERGE.out.merged
    )
    ch_versions = ch_versions.mix(DADA2_SEQUENCETABLE.out.versions)

    // Plot sequence distribution
    DADA2_SEQUENCEDISTRIBUTION (
        ch_ids,
        DADA2_SEQUENCETABLE.out.seq_table
    )
    ch_versions = ch_versions.mix(DADA2_SEQUENCEDISTRIBUTION.out.versions)

    // Remove Bimeras
    DADA2_REMOVEBIMERADENOVO (
        ch_ids,
        DADA2_SEQUENCETABLE.out.seq_table
    )
    ch_versions = ch_versions.mix(DADA2_REMOVEBIMERADENOVO.out.versions)

    // Track the read counts through the DADA2 stages
    DADA2_TRACKREADS (
        ch_prefix,
        ch_single_end,
        DADA2_SAMPLEINFERENCE.out.sample_inference,
        DADA2_REMOVEBIMERADENOVO.out.seq_table_nochim,
        DADA2_MERGE.out.merged,
        ch_filt_trim_csv
    )
    ch_versions = ch_versions.mix(DADA2_TRACKREADS.out.versions)

    // Create final outputs
    DADA2_FINALOUTPUTS (
        ch_prefix,
        ch_ids,
        DADA2_REMOVEBIMERADENOVO.out.seq_table_nochim
    )
    ch_versions = ch_versions.mix(DADA2_FINALOUTPUTS.out.versions)

    emit:
    fasta            = DADA2_FINALOUTPUTS.out.asv_fasta        // channel: [ val(prefix), asv.fa ]
    table            = DADA2_FINALOUTPUTS.out.asv_tsv          // channel: [ val(prefix), asv_final_table.tsv ]
    lca_input_table  = DADA2_FINALOUTPUTS.out.lca_input_tsv    // channel: [ val(prefix), lca_input.tsv ]
    quality_raw      = DADA2_PLOTQUALITYPROFILERAW.out.png     // channel: [ val(meta), ${meta.id}_raw*.png ]
    quality_filt     = DADA2_PLOTQUALITYPROFILETRIMMED.out.png // channel: [ val(meta), ${meta.id}_filtered*.png ]
    errors_plot      = DADA2_PLOTERRORS.out.png                // channel: [ errors_plot_*.png ]
    seq_distribution = DADA2_SEQUENCEDISTRIBUTION.out.png      // channel: [ ASV_seq_distribution.png ]
    track_reads      = DADA2_TRACKREADS.out.png                // channel: [ asv_track_reads.png ]
    versions         = ch_versions                             // channel: [ versions.yml ]
}
