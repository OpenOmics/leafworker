<div align="center">

  <h1 style="font-size: 250%">leafworker 🔬</h1>

  <b><i>An awesome snakemake pipeline to run leafcutter</i></b><br> 
  <a href="https://doi.org/10.5281/zenodo.15170953">
      <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.15170953.svg" alt="DOI">
  </a>
  <a href="https://github.com/OpenOmics/leafworker/releases">
    <img alt="GitHub release" src="https://img.shields.io/github/v/release/OpenOmics/leafworker?color=blue&include_prereleases">
  </a>
  <a href="https://hub.docker.com/repository/docker/skchronicles/leafcutter">
    <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/skchronicles/leafcutter">
  </a><br>
  <a href="https://github.com/OpenOmics/leafworker/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/leafworker/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/leafworker/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/leafworker/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/leafworker/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/leafworker?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/leafworker/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/leafworker">
  </a>

  <p>
    This is the home of the pipeline, leafworker. Its long-term goals: to make running leafcutter easier, more reproducible, and more scalable.
  </p>

</div>  


## Overview

Welcome to leafworker's documentation! This guide is the main source of documentation for users that are getting started with the [leafworker pipeline](https://github.com/OpenOmics/leafworker/). This pipeline is a wrapper around the [leafcutter<sup>0</sup>](https://github.com/davidaknowles/leafcutter/tree/master), a tool for the analysis of RNA-seq data to identify and quantify alternative splicing events. The pipeline is designed to be highly scalable, reproducible, and easy to use. If you use this pipeline, _please do not forget to cite [leafcutter](https://www.nature.com/articles/s41588-017-0004-9)_.

The **`./leafworker`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">leafworker <b>run</b></code>](usage/run.md)   
    Run the leafworker pipeline with your input files.

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">leafworker <b>unlock</b></code>](usage/unlock.md)  
    Unlocks a previous runs output directory.

</section>

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">leafworker <b>install</b></code>](usage/install.md)  
    Download remote reference files locally.


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">leafworker <b>cache</b></code>](usage/cache.md)  
    Cache remote software containers locally.  

</section>

**leafworker** is a pipeline to make running leafcutter easier, more reproducible, and more scalable. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of BAM files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/leafworker/issues).

## Contribute 

This site is a living document, created for and by members like you. leafworker is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/leafworker).

## Citation

If you use this software, please cite it as below:  

=== "BibTex"

    ```
    @software{Kuhn_OpenOmics_leafworker_2025,
      author       = {Kuhn, Skyler},
      title        = {OpenOmics/leafworker},
      month        = apr,
      year         = 2025,
      publisher    = {Zenodo},
      version      = {v0.1.0},
      doi          = {10.5281/zenodo.15170953},
      url          = {https://doi.org/10.5281/zenodo.15170953}
    }
    ```

=== "APA"

    ```
    Kuhn, S. (2025). OpenOmics/leafworker: v0.1.0. Zenodo. https://doi.org/10.5281/zenodo.15170953
    ```

For more citation style options, please visit the pipeline's [Zenodo page](https://doi.org/10.5281/zenodo.15170953).


## References

<sup>**0.** Li, Y. I., Knowles, D. A., Humphrey, J., Barbeira, A. N., Dickinson, S. P., Im, H. K., & Pritchard, J. K. (2018). Annotation-free quantification of RNA splicing using LeafCutter. Nature genetics, 50(1), 151–158. https://doi.org/10.1038/s41588-017-0004-9</sup>  
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
