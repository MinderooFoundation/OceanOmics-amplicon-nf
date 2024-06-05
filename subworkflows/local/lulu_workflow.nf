/*
 * Filter ZOTUs/ASVs with LULU
 */

include { BLAST_MATCHLIST   } from '../../modules/local/blast/matchlist/main'
include { BLAST_MAKEBLASTDB } from '../../modules/local/blast/makeblastdb/main'
include { LULU              } from '../../modules/local/lulu/main.nf'
include { GETLULUFASTA      } from '../../modules/local/custom/getlulufasta/main.nf'

workflow LULU_WORKFLOW {
    take:
    ch_fasta // channel: [ val(prefix), fasta ]
    ch_table // channel: [ val(prefix), table ]

    main:
    ch_versions = Channel.empty()

    // Make blast data base using fasta file
    BLAST_MAKEBLASTDB (
        ch_fasta
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    // Create blast list with blastn
    BLAST_MATCHLIST (
        ch_fasta,
        BLAST_MAKEBLASTDB.out.db
    )
    ch_versions = ch_versions.mix(BLAST_MATCHLIST.out.versions)

    // Run LULU
    LULU (
        ch_table,
        BLAST_MATCHLIST.out.txt
    )
    ch_versions = ch_versions.mix(LULU.out.versions)

    // Get a curated fasta file
    GETLULUFASTA (
        ch_fasta,
        LULU.out.curated_table
    )

    emit:
    curated_table = LULU.out.curated_table // channel: [ val(prefix), curated_table.tab ]
    curated_fasta = GETLULUFASTA.out.fasta // channel: [ val(prefix), curated.fasta ]
    versions      = ch_versions            // channel: [ versions.yml ]
}
