# MinderooFoundation/OceanOmics-amplicon-nf: Parameters

## Introduction

All the parameters in the pipeline can be set in a config file, or they can be set at the command line.

### General parameters

- [input](#input) - Input sample sheet .csv file (mandatory columns if skipping demultiplexing: sample, fastq_1, fastq_2)
  (mandatory columns if demultiplexing: sample, fastq_1, fastq_2, fw_index, rv_index, fw_primer, rv_primer, fw_no, rv_no)
- [bind_dir](#bind_dir) - Directory to bind to the Docker images (e.g., /scratch)
- [dbfiles](#dbfiles) - List of Blast database files in quote marks (e.g., "path/to/db/\*")
- [filter_table](#filter_table) - Optional .csv file to filter out unwanted taxa. First column should be called `level`, and second column should be called `name`. For example, one row might have domain in the `level` column and Bacteria in the `name` column.

### Demultiplex parameters

- [raw_data](#raw_data) - Raw data file/s in quote marks; use {} if your data is paired end (e.g., "path/to/fqs/prefix*{R1,R2}*.fq.gz")
- [cutadapt_error](#cutadapt_error) - default = 0.15
- [cutadapt_min_len](#cutadapt_min_len) - default = 1
- [obitools3_min_len](#obitools3_min_len) - default = 8
- [obi3_demux](#obi3_demux) - Demultiplex with obitools3 instead of Cutadapt; default = false
- [ulimit](#ulimit) - Increase this value if you run into a 'Too many open files error' during Cutadapt; default = 10000

### ZOTU parameters

- [min_quality](#min_quality) - default = 10
- [min_align_leng](#min_align_leng) - default = 6
- [min_size](#min_size) - default = 8
- [id](#id) - default = 0.97
- [mode](#mode) - default = "vsearch", can also be "usearch32", or "usearch64"
- [usearch64](#usearch64) - Full path to where usearch64 bit version is stored locally

### ASV parameters

- [pooled](#pooled) - default = "true"; can also be "false" or "pseudo"
- [min_overlap](#min_overlap) - default = 12
- [max_mismatch](#max_mismatch) - default = 0

### Blast parameters

- [blast_task](#blast_task) - default = "blastn"
- [perc_identity](#perc_identity) - default = 95
- [evalue](#evalue) - default = "1e-3"
- [best_hit_score_edge](#best_hit_score_edge) - default = 0.05
- [best_hit_overhang](#best_hit_overhang) - default = 0.25
- [qcov](#qcov) - default = 100
- [max_tar_seq](#max_tar_seq) - default = 10

### LULU parameters

- [min_match_lulu](#min_match_lulu) - default = 84

### LCA parameters

- [lca_qcov](#lca_qcov) - default = 100
- [lca_pid](#lca_pid) - default = 97
- [lca_diff](#lca_diff) - default = 1

### Parameters to skip steps

- [skip_demux](#skip_demux) - default = false
- [skip_asvs](#skip_asvs) - default = false
- [skip_zotus](#skip_zotus) - default = false
- [skip_lulu](#skip_lulu) - default = false
