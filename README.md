# Have Freyja do all the things

`alltheFreyjas` is a simple Nextflow wrapper around commands from the [`freyja` suite of commands](https://github.com/andersen-lab/Freyja), which quantify the SARS-CoV-2 lineages in wastewater samples. It uses the [`staphb`](https://hub.docker.com/r/staphb/freyja) docker image for `freyja` to ensure reproducibility, and spreads `freyja` executions across available cores.

### Quick start

To download the pipeline, simply git clone it into an empty directory of your choice with:

```
git clone https://github.com/nrminor/alltheFreyjas.git .
```

Assuming you have [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/) installed, run it on a directory of SARS-CoV-2 BAM files with:

```
nextflow run main.nf --bam_dir BAM_DIR
```

The pipeline also supports Singularity/Apptainer so that it can be run on HPC clusters.
