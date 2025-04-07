<div align="center">
   
  <h1>leafworker 🔬</h1>
  
  **_An awesome snakemake pipeline to run leafcutter_**

  [![tests](https://github.com/OpenOmics/leafworker/workflows/tests/badge.svg)](https://github.com/OpenOmics/leafworker/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/leafworker/workflows/docs/badge.svg)](https://github.com/OpenOmics/leafworker/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/leafworker?color=brightgreen)](https://github.com/OpenOmics/leafworker/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/leafworker)](https://github.com/OpenOmics/leafworker/blob/main/LICENSE) 
  
  <i>
    This is the home of the pipeline, leafworker. Its long-term goals: to make running leafcutter easier, more reproducible, and more scalable.
  </i>
</div>

## Overview

Welcome to leafworker! Before getting started, we highly recommend reading through [leafworker's documentation](https://openomics.github.io/leafworker/). This pipeline is a wrapper around the [leafcutter<sup>0</sup>](https://github.com/davidaknowles/leafcutter/tree/master), a tool for the analysis of RNA-seq data to identify and quantify alternative splicing events. The pipeline is designed to be highly scalable, reproducible, and easy to use. If you use this pipeline, _please do not forget to cite [leafcutter](https://www.nature.com/articles/s41588-017-0004-9)_.

The **`./leafworker`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>leafworker <b>run</b></code>](https://openomics.github.io/leafworker/usage/run/): Run the leafworker pipeline with your input files.
 * [<code>leafworker <b>unlock</b></code>](https://openomics.github.io/leafworker/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>leafworker <b>install</b></code>](https://openomics.github.io/leafworker/usage/install/): Download reference files locally.
 * [<code>leafworker <b>cache</b></code>](https://openomics.github.io/leafworker/usage/cache/): Cache software containers locally.

**leafworker** is a pipeline to make running leafcutter easier, more reproducible, and more scalable. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of BAM files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/leafworker/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/leafworker/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/leafworker/issues).

## Dependencies

**Requires:** `singularity>=3.5`  `snakemake<=7.32.4`

At the current moment, the pipeline only has two dependencies: snakemake and singularity. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline relies on versioned images from [DockerHub](https://hub.docker.com/repository/docker/skchronicles/leafcutter). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation

Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/leafworker.git
# Change your working directory
cd leafworker/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./leafworker -h
```

## Contribute 

This site is a living document, created for and by members like you. leafworker is maintained by the members of OpenOmics and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/leafworker).

## Cite

If you use this software, please cite it as below:  

<details>
  <summary><b><i>@BibText</i></b></summary>
 
```text
Citation coming soon!
```

</details>

<details>
  <summary><b><i>@APA</i></b></summary>

```text
Citation coming soon!
```

</details>


## References

<sup>**0.** Li, Y. I., Knowles, D. A., Humphrey, J., Barbeira, A. N., Dickinson, S. P., Im, H. K., & Pritchard, J. K. (2018). Annotation-free quantification of RNA splicing using LeafCutter. Nature genetics, 50(1), 151–158. https://doi.org/10.1038/s41588-017-0004-9</sup>   
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
