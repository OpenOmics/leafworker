# Functions and rules for running IsoformSwitchAnalyzer
# Standard Library
import textwrap
# Local imports
from scripts.common import (
    allocated
)


# Data processing ISA rules
rule isoformswitchanalyzer_bam2fastq:
    """
    Data-processing step to convert BAM files into FastQ files.
    To get transcript counts salmon is run. Salmon only takes
    BAM files as input if they were aligned against the
    transcriptome. The pipeline takes genomic BAM files as
    input, so we will convert the BAM files into paired-end
    FastQ files and run Salmon with the paired-end FastQ
    files as input.
    @Input:
        Input BAM file (scatter-per-sample)
    @Output:
        Paired-end FastQ files
    """
    input:
        bam   = join(workpath, "{name}.bam"),
    output:
        r1 = join(workpath, "fastqs", "{name}.R1.fastq.gz"),
        r2 = join(workpath, "fastqs", "{name}.R2.fastq.gz"),
    params:
        rname = "bam2fqs",
    resources:
        mem   = allocated("mem",  "isoformswitchanalyzer_bam2fastq", cluster),
        time  = allocated("time", "isoformswitchanalyzer_bam2fastq", cluster),
    threads: int(allocated("threads", "isoformswitchanalyzer_bam2fastq", cluster))
    container: config["images"]["isoformswitchanalyzer"]
    shell: """
    # Creates paired-end FastQ files while 
    # discarding singletons, supplementary, 
    # and secondary reads using samtools
    samtools collate -@ {threads} -u -O {input.bam} \\
    | samtools fastq \\
        -1 {output.r1} \\
        -2 {output.r2} \\
        -0 /dev/null \\
        -s /dev/null \\
        -n
    """
