# MinderooFoundation/OceanOmics-amplicon-nf: Parameters

## Introduction

All the parameters in the pipeline can be set in a config file, or they can be set at the command line.

### General parameters

- [input](#input) - Input sample sheet .csv file (mandatory columns if skipping demultiplexing: sample, fastq_1, fastq_2)
  (mandatory columns if demultiplexing: sample, fastq_1, fastq_2, fw_index, rv_index, fw_primer, rv_primer, fw_no, rv_no)
- [outdir](#outdir) - Directory where output files will be published
- [binddir](#bind_dir) - Directory to bind to the Docker images (e.g., /scratch)
- [dbfiles](#dbfiles) - List of Blast database files in quote marks (e.g., "path/to/db/\*")
- [filter_table](#filter_table) - Optional .csv file to filter out unwanted taxa. First column should be called `level`, and second column should be called `name`. For example, one row might have domain in the `level` column and Bacteria in the `name` column.
- [fw_primer](#fw_primer) - Your forward primer. Can be a semi-colon seperated list of primers.
- [rv_primer](#rv_primer) - Your reverse primer. Can be a semi-colon seperated list of primers.

### Parameters to automate choosing other parameters

- [assay](#assay) - The assay of your data if you want certain parameters chosen automatically. Currently supports `16SFish`, `16SMam`, `MiFish`, `12SV5`, and `COILeray`. The `-c` option can be used to provide a custom map with other assays. More information can be found [here](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf/blob/master/docs/custom_config.md)

### Demultiplex parameters

- [raw_data](#raw_data) - Raw data file/s in quote marks; use {} if your data is paired end (e.g., "path/to/fqs/prefix*{R1,R2}*.fq.gz")
- [cutadapt_error](#cutadapt_error) - default = 0
- [cutadapt_min_len](#cutadapt_min_len) - default = 1
- [obi3_min_len](#obi3_min_len) - default = 8
- [obi3_demux](#obi3_demux) - Demultiplex with obitools3 instead of Cutadapt; default = false
- [ulimit](#ulimit) - Increase this value if you run into a 'Too many open files error' during Cutadapt; default = 40000
- [demux_udi](#demux_udi) - Demultiplex unique dual-indexes (The default method assumes unique combinatorial indexes)

### FAIRe parameters

- [faire_mode](#faire_mode) - Turns on FAIRe mode that and adds taxa results to the metadata while following FAIRe standards (https://fair-edna.github.io/index.html)
- [faire_metadata](#faire_metadata) - Metadata that needs to follow FAIRe standards. It's recommended to use FAIRe-ator https://github.com/FAIR-eDNA/FAIRe-ator/tree/main

### Trim parameters

- [primertrim_error](#primertrim_error) - The error rate to allow mismatches when trimming primer sequences, default = 0
- [seqtk_trim](#seqtk_trim) - Option to use seqtk to trim BPs after a length cutoff
- [seqtk_length](#seqth_length) - Cutoff length to used with seqtk_trim. Every BP after the cutoff will be trimmed off; default = 180
- [fastp_trim](#fastp_trim) - Option to use fastp to trim BPs from each end using a sliding window
- [fastp_front_window_size](#fastp_front_window_size) - sliding window size when trimming from the left side of a read; default = 4
- [fastp_front_mean_quality](#fastp_front_mean_quality) - mean quality cutoff when trimming from the left side of a read; default = 20
- [fastp_tail_window_size](#fastp_tail_window_size) - sliding window size when trimming from the right side of a read; default = 4
- [fastp_tail_mean_quality](#fastp_tail_mean_quality) - mean quality cutoff when trimming from the right side of a read; default = 20

### ZOTU parameters

- [zotu_min_quality](#zotu_min_quality) - default = 10
- [zotu_min_align_leng](#zotu_min_align_leng) - default = 6
- [zotu_min_size](#zotu_min_size) - default = 8
- [zotu_id](#zotu_id) - default = 0.97
- [zotu_mode](#zotu_mode) - default = "vsearch", can also be "usearch32", or "usearch64"
- [usearch64](#usearch64) - Full path to where usearch64 bit version is stored locally

### ASV parameters

- [asv_pooled](#asv_pooled) - default = "true"; can also be "false" or "pseudo"
- [asv_min_overlap](#asv_min_overlap) - default = 12
- [asv_max_mismatch](#asv_max_mismatch) - default = 0

### Blast parameters

- [blast_task](#blast_task) - default = "blastn"
- [blast_pid](#blast_pid)
- [blast_evalue](#blast_evalue) - default = "0.001"
- [blast_best_hit_score_edge](#blast_best_hit_score_edge)
- [blast_best_hit_overhang](#blast_best_hit_overhang)
- [blast_qcov](#blast_qcov) - default = 100
- [blast_max_tar_seq](#blast_max_tar_seq) - default = 999999

### LULU parameters

- [lulu_min_match](#lulu_min_match) - default = 84

### LCA parameters

- [lca_qcov](#lca_qcov) - default = 100
- [lca_pid](#lca_pid) - default = 90
- [lca_diff](#lca_diff) - default = 1
- [lca_with_fishbase](#lca_with_fishbase). Use alternate LCA script that checks Fishbase, then Worms, then NCBI for lineage information. - default = false

### Parameters to skip steps

- [skip_demux](#skip_demux). Skip demultiplexing. - default = false
- [skip_primer_trim](#skip_primer_trim). Skip primer trimming from the 5' end of reads. This step discards reads that are missing a primer - default = false
- [skip_primer_leftover_trim](#skip_primer_leftover_trim). By default, the pipeline trims reverse complemented primers from the 3' end of reads while also trimming off any primer-dimers. This parameter skips that step. - default = false
- [skip_asvs](#skip_asvs). This will prevent the pipeline from creating ASVs. - default = false
- [skip_zotus](#skip_zotus). This will prevent the pipeline from creating ZOTUs. - default = false
- [skip_lulu](#skip_lulu). This will skip LULU. - default = false
- [skip_lulu_comparison](#skip_lulu_comparison). By default, the pipeline produces phyloseq objects with and without LULU curation. If you're not skipping LULU, this will prevent creating non-LULU phyloseq objects. - default = false
- [skip_classification](#skip_classification). This will skip the Blastn/LCA part of the pipeline. - default = false
- [skip_filter](#skip_proportionalfilter). This will skip the proportional filter part of the pipeline. Proportional filter is a filter step using the filter method described here (https://essopenarchive.org/doi/full/10.22541/au.169956117.76591919) - default = false
- [normalise_pid](#normalise_pid). The LCA_withFishbase can normalise the percent ID by multiplying that value with the query coverage. This option turn that feature on - default = false
- [use_bitwise](#use_bitwise). By default, the LCA_withFishbase script uses percent identity when calculating LCA. This option tells the script to use bitwise instead of precent identity. Default = false
- [start_from_blast](#start_from_blast). This will start the pipeline at the BLAST step. If using this option, you must also provide the "fasta" and "otu_table". Default = false
- [start_from_lca](#start_from_lca). This will start the pipeline at the LCA step. If using this option, you must also provide the "fasta", "otu_table", and "blast_results". Default = false
- [fasta](#fasta). ASV/ZOTU fasta file to be used when starting from the BLAST or LCA steps.
- [otu_table](#otu_table). ASV/ZOTU count table to be used when starting from the BLAST or LCA steps.
- [blast_results](#blast_results). BLAST results file to be used when starting at the LCA step. BLAST should have been run with "-outfmt '6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp'"
