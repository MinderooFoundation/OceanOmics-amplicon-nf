## Introduction

This pipeline is used to create ASVs and ZOTUs from eDNA amplicon data, assign taxonomy to those ASVs/ZOTUs and finally produce phyloseq objects.

**OceanOmics-amplicon-nf** creates a phyloseq object from eDNA amplicon data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
3. Demultiplex with Cutadapt ([`Cutadapt`])(https://cutadapt.readthedocs.io/en/stable/)
4. Optionally demultiplex with Obitools3 ([`Obitools3`])(https://git.metabarcoding.org/obitools/obitools3)
5. Additional QC with Seqkit Stats ([`Seqkit`])(https://bioinf.shenwei.me/seqkit/usage/)
6. Create ASVs with DADA2 ([`DADA2`](https://www.bioconductor.org/packages/release/bioc/html/dada2.html))
7. Create ZOTUs with VSEARCH ([`VSEARCH`](https://manpages.ubuntu.com/manpages/bionic/man1/vsearch.1.html))
8. Optionally create ZOTUs with USEARCH ([`USEARCH`])(https://www.drive5.com/usearch/) 
9. Cutate ASVs/ZOTUs with LULU ([`LULU`](https://github.com/tobiasgf/lulu))
10. Assign taxonomy with blastn ([`blastn`](https://www.ncbi.nlm.nih.gov/books/NBK279691/))
11. Lowest Common Ancestor ([`LCA`](https://github.com/mahsa-mousavi/eDNAFlow/tree/master/LCA_taxonomyAssignment_scripts))
12. Phyloseq object creation ([`phyloseq`](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline with command:

   ```bash
   git clone git@github.com:MinderooFoundation/OceanOmics-amplicon-nf.git
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```bash
   nextflow run path/to/main.nf --input samplesheet.csv --outdir <OUDIR> --bind_dir <BINDDIR> --dbfiles "<BLASTDBFILES>" -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> 
   ```
## Documentation

The **OceanOmics-amplicon-nf** pipeline comes with documentation about the pipeline [usage](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf/blob/master/docs/usage.md), [parameters](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf/blob/master/docs/parameters.md) and [output](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf/blob/master/docs/usage.md).

## Credits

This pipeline incorporates aspects of eDNAFlow, which was written by Mahsa Mousavi. **OceanOmics-amplicon-nf** was written by Adam Bennett. Other people who have contributed to this pipeline include Sebastian Rauschert (conceptualisation), Philipp Bayer, and Jessica Pearce. This pipeline was built using the nf-core template.

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/ednaflow for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

