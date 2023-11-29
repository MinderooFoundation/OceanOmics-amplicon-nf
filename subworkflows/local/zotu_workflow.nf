/*
 * ZOTU Creation with VSEARCH, USEARCH32, or USEARCH64
 */

include { FASTQRELABEL            } from '../../modules/local/custom/fastqrelabel/main'
include { ADAPTERREMOVAL          } from '../../modules/local/adapterremoval/main'
include { FASTQCONCAT             } from '../../modules/local/custom/fastqconcat/main'
include { FASTQTOFASTA            } from '../../modules/local/custom/fastqtofasta/main'
include { USEARCH32_FASTXUNIQUES  } from '../../modules/local/usearch32/fastxuniques/main'
include { USEARCH32_UNOISE3       } from '../../modules/local/usearch32/unoise3/main'
include { USEARCH32_OTUTAB        } from '../../modules/local/usearch32/otutab/main'
include { USEARCH64_FASTXUNIQUES  } from '../../modules/local/usearch64/fastxuniques/main'
include { USEARCH64_UNOISE3       } from '../../modules/local/usearch64/unoise3/main'
include { USEARCH64_OTUTAB        } from '../../modules/local/usearch64/otutab/main'
include { VSEARCH_DEREPFULLLENGTH } from '../../modules/local/vsearch/derepfulllength/main'
include { VSEARCH_CLUSTERUNOISE   } from '../../modules/local/vsearch/clusterunoise/main'
include { VSEARCH_UCHIME3DENOVO   } from '../../modules/local/vsearch/uchime3denovo/main'
include { VSEARCH_USEARCHGLOBAL   } from '../../modules/local/vsearch/usearchglobal/main'

workflow ZOTU_WORKFLOW {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()
    ch_prefix   = "zotu"

    // Relabel fastq files
    FASTQRELABEL (
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQRELABEL.out.versions.first())

    // Quality filtering and merge if reads are paired-end
    ADAPTERREMOVAL (
        FASTQRELABEL.out.reads
    )
    ch_versions = ch_versions.mix(ADAPTERREMOVAL.out.versions.first())

    // We are collecting the samples together because the rest of the VSEARCH workflow needs the samples collected
    ch_reads_filt = ADAPTERREMOVAL.out.reads.map { it[1] }.collect()

    // Concatenate all the fastq files
    FASTQCONCAT (
        ch_reads_filt
    )
    ch_versions = ch_versions.mix(FASTQCONCAT.out.versions)

    // Convert fastq to fasta
    FASTQTOFASTA (
        FASTQCONCAT.out.reads
    )
    ch_versions = ch_versions.mix(FASTQTOFASTA.out.versions)

    if (params.zotu_mode == "vsearch") {
        // Dereplicate
        VSEARCH_DEREPFULLLENGTH (
            FASTQTOFASTA.out.reads
        )
        ch_versions = ch_versions.mix(VSEARCH_DEREPFULLLENGTH.out.versions)

        // Denoise
        VSEARCH_CLUSTERUNOISE (
            VSEARCH_DEREPFULLLENGTH.out.reads
        )
        ch_versions = ch_versions.mix(VSEARCH_CLUSTERUNOISE.out.versions)

        // Remove Chimeras
        VSEARCH_UCHIME3DENOVO (
            ch_prefix,
            VSEARCH_CLUSTERUNOISE.out.reads
        )
        ch_versions = ch_versions.mix(VSEARCH_UCHIME3DENOVO.out.versions)

        // Global pairwise alignment
        VSEARCH_USEARCHGLOBAL (
            FASTQTOFASTA.out.reads,
            VSEARCH_UCHIME3DENOVO.out.zotu_fasta
        )
        ch_versions = ch_versions.mix(VSEARCH_USEARCHGLOBAL.out.versions)

        ch_derep_fasta   = VSEARCH_DEREPFULLLENGTH.out.reads
        ch_denoise_fasta = VSEARCH_UCHIME3DENOVO.out.zotu_fasta
        ch_zotu_fasta    = VSEARCH_UCHIME3DENOVO.out.zotu_fasta
        ch_zotu_tsv      = VSEARCH_USEARCHGLOBAL.out.zotu_tsv
        ch_lca_input     = VSEARCH_USEARCHGLOBAL.out.lca_input_tsv

    } else if (params.zotu_mode == "usearch32") {
        // Dereplicate
        USEARCH32_FASTXUNIQUES (
            FASTQTOFASTA.out.reads
        )
        ch_versions = ch_versions.mix(USEARCH32_FASTXUNIQUES.out.versions)

        // Denoise
        USEARCH32_UNOISE3 (
            ch_prefix,
            USEARCH32_FASTXUNIQUES.out.reads
        )
        ch_versions = ch_versions.mix(USEARCH32_UNOISE3.out.versions)

        // Global pairwise alignment
        USEARCH32_OTUTAB (
            FASTQTOFASTA.out.reads,
            USEARCH32_UNOISE3.out.zotu_fasta
        )
        ch_versions = ch_versions.mix(USEARCH32_OTUTAB.out.versions)

        ch_derep_fasta   = USEARCH32_FASTXUNIQUES.out.reads
        ch_denoise_fasta = USEARCH32_UNOISE3.out.zotu_fasta
        ch_zotu_fasta    = USEARCH32_UNOISE3.out.zotu_fasta
        ch_zotu_tsv      = USEARCH32_OTUTAB.out.zotu_tsv
        ch_lca_input     = USEARCH32_OTUTAB.out.lca_input_tsv

    } else if (params.zotu_mode == "usearch64") {
        usearch64 = params.usearch64

        // Dereplicate
        USEARCH64_FASTXUNIQUES (
            usearch64,
            FASTQTOFASTA.out.reads
        )
        ch_versions = ch_versions.mix(USEARCH64_FASTXUNIQUES.out.versions)

        // Denoise
        USEARCH64_UNOISE3 (
            usearch64,
            USEARCH64_FASTXUNIQUES.out.reads
        )
        ch_versions = ch_versions.mix(USEARCH64_UNOISE3.out.versions)

        // Global pairwise alignment
        USEARCH64_OTUTAB (
            usearch64,
            FASTQTOFASTA.out.reads,
            USEARCH64_UNOISE3.out.zotu_fasta
        )
        ch_versions = ch_versions.mix(USEARCH64_OTUTAB.out.versions)

        ch_derep_fasta   = USEARCH64_FASTXUNIQUES.out.reads
        ch_denoise_fasta = USEARCH64_UNOISE3.out.zotu_fasta
        ch_zotu_fasta    = USEARCH64_UNOISE3.out.zotu_fasta
        ch_zotu_tsv      = USEARCH64_OTUTAB.out.zotu_tsv
        ch_lca_input     = USEARCH64_OTUTAB.out.lca_input_tsv
    }

    emit:
    fasta           = ch_zotu_fasta          // channel: [ val(prefix), zotu.fa ]
    table           = ch_zotu_tsv            // channel: [ val(prefix), zotu_table.tsv ]
    lca_input_table = ch_lca_input           // channel: [ val(prefix), lca_input.tsv ]
    versions        = ch_versions            // channel: [ versions.yml ]
}
