
## Introduction

YAPGAP is a fork of **nf-core/genomeannotator** (a bioinformatics best-practice analysis pipeline for the annotation of eukaryote genomes.

YAPGAP was developed some challenges with the speed and issues with genomeannotator, so some tools have been replaced, however, overall it is copy of the **nf-core/genomeannotator** and we gratefull for their pipeline and therefore some to the content here remains the  sames the original pipeline 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

1. Preprocess assembly (filter out small contigs, clean sequence names)
2. Align evidences (proteins, transcripts, RNAseq)
3. Convert alignments to annotation support
4. Build gene models using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) and optionally [PASA](https://github.com/PASApipeline/PASApipeline)
5. Compute consensus gene build using [EvidenceModeler](https://evidencemodeler.github.io/)

Optional steps include de-novo transcriptome assembly ([Trinity](https://github.com/trinityrnaseq/trinityrnaseq)) and annotation mapping from related genomes ([Satsuma2](https://github.com/bioinfologics/satsuma2) and [Kraken](https://github.com/GrabherrGroup/kraken)).

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/genomeannotator -profile test,YOURPROFILE --outdir <OUTDIR>
    ```

    Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

    > * The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    ```console
    nextflow run nf-core/genomeannotator --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

## Dependencies

+ [Genome Assembly Annotation Service (GAAS)](https://github.com/NBISweden/GAAS)
  - Assembly preprocessing
    + filtering sequences by size


## Credits

nf-core/genomeannotator was originally written by Marc P. Hoeppner
