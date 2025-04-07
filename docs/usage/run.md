# <code>leafworker <b>run</b></code>

## 1. About 

The `leafworker` executable is composed of several inter-related sub commands. Please see `leafworker -h` for all available options.

This part of the documentation describes options and concepts for <code>leafworker <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running leafworker pipeline. 

Setting up the leafworker pipeline is fast and easy! In its most basic form, <code>leafworker <b>run</b></code> only has *three required inputs*.

## 2. Synopsis

```text
$ leafworker run [--help] \
      [--dry-run] [--job-name JOB_NAME] [--mode {{slurm,local}}] \
      [--sif-cache SIF_CACHE] [--singularity-cache SINGULARITY_CACHE] \
      [--silent] [--threads THREADS] [--tmp-dir TMP_DIR] \
      [--batch-id BATCH_ID] [--groups GROUPS] [--contrasts CONTRASTS] \
      --input INPUT [INPUT ...] \
      --output OUTPUT \
      --gtf GTF
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of BAM (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, and a GTF file to define the transcriptome via `--gtf`.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input BAM file(s).**  
> *type: BAM file(s)*  
> 
> Input BAM files to process. BAM files for one or more samples can be provided. Multiple input BAM files should be seperated by a space. Globbing for multiple file is also supported! This makes selecting BAM files easy.
> 
> ***Example:*** `--input .tests/*.bam`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/leafworker_out`

---  
  `--gtf GTF`
> **Annotation file in GTF format.**   
> *type: GTF file*
>   
> This annotation file is used by leafcutter to identify splice junctions within the transcriptome. The input annotation file can be compressed with gzip.
> 
> ***Example:*** `--gtf .tests/gencode.v19.annotation.gtf.gz`

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. 

  `--batch-id BATCH_ID`
> **Unique identifer to associate with a batch of samples.**   
> *type: string*  
> *default: None*   
>
> This option can be provided to ensure that differential splicing output files are not over-written between runs of the pipeline after updating the group file with additional covariates or dropping samples. By default, project-level files in "differential_splicing" could get over-written between pipeline runs if this option is not provided. The output directory name for a given contrast will resolve to `{group1}_vs_{group2}` within the **differential_splicing** folder. As so, if the groups file is updated to remove samples or add additional covariates without updating the group names, it could over write the previous analyses output files. Any identifer provided to this option will be used to create a sub directory in the "differential_splicing" folder. This ensures project-level files (which are unique) will not get over written. With that being said, it is always a good idea to provide this option. A unique batch id should be provided between runs. This batch id should be composed of alphanumeric characters and it should not contain a white space or tab characters. Here is a list of valid or acceptable characters: `aA-Zz`, `0-9`, `-`, `_`.
> 
> ***Example:*** `--batch-id 2025_04_08`

---  
  `--groups GROUPS`
> **Groups file containing sample metadata.**   
> *type: TSV file*
>   
> This tab delimited file is used to pair each sample to a group. Group information is used to setup comparsions across different groups of samples to find differential splicing. This tab-delimited file consists of two columns containing the base name of each sample and the name of their group(s), where  multiple groups seperated can be seperated by a comma. The header of this file needs to be `Sample` for the column containing the sample base names and `Group` for the column with the group names. The base name of a given sample can be determined by removing its file extension from the sample's bam file.  For example: *WT_S4.bam* becomes *WT_S4* in the groups file. A group can represent anything from a timepoint, an experimental condition, a treatment, or disease state, etc.  Please note that groups are should be composed of alphanumeric characters, cannot startwith a number, and cannot contain any `-` characters. Groups can contain `_` characters. 
> 
> In the groups file the **1st & 2nd columns are required**, where the **1st** column must be called `Sample` and the **2nd** column must be called `Group`. Any additional columns after the `Sample` & `Group` columns are optional. Any extra columns will be used for controlling covariates. This could be batches, sex, age, etc. *Please note that any numeric values will be modeled as continuous variables*. As so, please ensure any categorical values start with a letter. To perform differential splicing  analyses, a groups & contrast file must be provided.
>
> *Here are the contents of example groups file:*
> ```
> Sample	Group	Sex
> WT_S1	G1,G3	M
> WT_S2	G1,G3	F
> WT_S3	G1,G3	M
> WT_S4	G1,G3	F
> WT_S5	G1,G4	M
> WT_S6	G2,G4	F
> WT_S7	G2,G4	M
> WT_S8	G2,G4	F
> WT_S9	G2,G4	M
> WT_S10	G2,G4	F
> ```
> 
> ***Example:*** `--groups .tests/groups.tsv`

---  
  `--contrasts CONTRASTS`
> **Contrasts file containing comparisons to make.**   
> *type: TSV file*
>   
> This tab delimited file is used to setup comparisons within different groups of samples. Please see the `--groups` option above for more information about how to define groups within a set of samples. The contrasts file consists of two columns containing the names of each group to compare. *Please note this file has no column names or header.* The group names defined in this file must also exist in the groups file. To perform differential splicing analysis, a groups and contrast file must be provided together.
>
> *Here are the Contents of example contrasts file:*
> ```
> G2	G1
> G4	G3
> ```
> 
> ***Example:*** `--contrasts .tests/contrasts.tsv`


### 2.3 Orchestration options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*   
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local. 
> 
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running leafworker in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode. 
> 
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*  
> *default: pl:leafworker*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:leafworker".
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `leafworker cache` subcommand can be used to create a local SIF cache. Please see `leafworker cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running leafworker with this option when ever possible.
> 
> ***Example:*** `--sif-cache /data/$USER/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`


---  
  `--tmp-dir TMP_DIR`   
> **Max number of threads for each process.**  
> *type: path*  
> *default: None*
> 
> Path on the file system for writing temporary output files. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in `/data/scratch`. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
> 
> ***Example:*** `--tmp-dir /data/scratch/$USER`

### 2.4 Miscellaneous options 

Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example

```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load snakemake/7.22.0-ufanewz

# Step 2A.) Dry-run the pipeline
./leafworker run --input .tests/*.bam \
                  --output /data/$USER/output \
                  --gtf .tests/gencode.v19.annotation.gtf.gz \
                  --mode slurm \
                  --dry-run

# Step 2B.) Run the leafworker pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode.
./leafworker run --input .tests/*.bam \
                  --output /data/$USER/output \
                  --gtf .tests/gencode.v19.annotation.gtf.gz \
                  --mode slurm
```