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
- [fw_primer](#fw_primer) - Your forward primer. Reads will be checked to see if they're reverse complemented.
- [rv_primer](#rv_primer) - Your reverse primer. Reads will be checked to see if they're reverse complemented.

### Parameters to automate choosing other parameters

- [assay](#assay) - The assay of your data if you want certain parameters chosen automatically. Currently supports `16SFish`, `16SMam`, `MiFish`, `12SV5`, and `COILeray`. The `-c` option can be used to provide a custom map with other assays. More information can be found [here](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf/blob/master/docs/custom_config.md)

### Demultiplex parameters

- [raw_data](#raw_data) - Raw data file/s in quote marks; use {} if your data is paired end (e.g., "path/to/fqs/prefix*{R1,R2}*.fq.gz")
- [cutadapt_error](#cutadapt_error) - default = 0
- [cutadapt_min_len](#cutadapt_min_len) - default = 1
- [obi3_min_len](#obi3_min_len) - default = 8
- [obi3_demux](#obi3_demux) - Demultiplex with obitools3 instead of Cutadapt; default = false
- [ulimit](#ulimit) - Increase this value if you run into a 'Too many open files error' during Cutadapt; default = 40000

### Trim parameters

- [primertrim_error](#primertrim_error) - The error rate to allow mismatches when trimming primer sequences, default = 2
- [seqtk_trim](#seqtk_trim) - Option to use seqtk to trim BPs after a length cutoff
- [seqtk_length](#seqth_length) - Cutoff length to used with seqtk_trim. Every BP after the cutoff will be trimmed off; default = 180

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

### Parameters to skip steps

- [skip_demux](#skip_demux). Skip demultiplexing. - default = false
- [skip_primertrim](#skip_primertrim). Skip primer trimming. - default = false
- [skip_3end_trim](#skip_3end_trim). By default, the pipeline trims reverse complemented primers from the 3' end of reads, this skips that step. - default = false
- [skip_asvs](#skip_asvs). This will prevent the pipeline from creating ASVs. - default = false
- [skip_zotus](#skip_zotus). This will prevent the pipeline from creating ZOTUs. - default = false
- [skip_lulu](#skip_lulu). This will skip LULU. - default = false
- [skip_lulu_comparison](#skip_lulu_comparison). By default, the pipeline produces phyloseq objects with and without LULU curation. If you're not skipping LULU, this will prevent creating non-LULU phyloseq objects. - default = false
- [skip_classification](#skip_classification). This will skip the Blastn/LCA part of the pipeline. - default = false
- [skip_nesterfilter](#skip_nesterfilter). This will skip the Nesterfilter part of the pipeline. Nesterfilter is a filter step using the filter method described here (https://essopenarchive.org/doi/full/10.22541/au.169956117.76591919) - default = false
