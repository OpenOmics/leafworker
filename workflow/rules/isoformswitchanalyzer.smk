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
        r1 = temp(join(workpath, "fastqs", "{name}.R1.fastq.gz")),
        r2 = temp(join(workpath, "fastqs", "{name}.R2.fastq.gz")),
    params:
        rname  = "bam2fqs",
        tmpdir = join(workpath, "temp"),
    resources:
        mem   = allocated("mem",  "isoformswitchanalyzer_bam2fastq", cluster),
        time  = allocated("time", "isoformswitchanalyzer_bam2fastq", cluster),
    threads: int(allocated("threads", "isoformswitchanalyzer_bam2fastq", cluster))
    container: config["images"]["isoformswitchanalyzer"]
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    export TMPDIR="${{tmp}}"

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


rule isoformswitchanalyzer_salmon_index:
    """
    Reference building step to create a index of the
    transcriptome for salmon. Salmon requires a pre-built 
    index of the transcriptome to perform quantification.
    @Input:
        Transcriptomic FASTA file (singleton)
    @Output:
        Salmon Index
    """
    input:
        transcripts = quantify_transcripts,
    output:
        index = join(workpath, "temp", "salmon_index", "seq.bin"),
    params:
        rname  = "salmonidx",
        prefix = join(workpath, "temp", "salmon_index"),
        tmpdir = join(workpath, "temp"),
    resources:
        mem   = allocated("mem",  "isoformswitchanalyzer_salmon_index", cluster),
        time  = allocated("time", "isoformswitchanalyzer_salmon_index", cluster),
    threads: int(allocated("threads", "isoformswitchanalyzer_salmon_index", cluster))
    container: config["images"]["isoformswitchanalyzer"]
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    export TMPDIR="${{tmp}}"

    # Build an index for running salmon
    salmon index \\
        --threads {threads} \\
        -t {input.transcripts} \\
        -i {params.prefix} \\
        -k 31 \\
        --gencode
    """


rule isoformswitchanalyzer_salmon_quant:
    """
    Data-processing step to get counts using Salmon. Salmon
    is an extremely fast and bias-aware method for transcript
    quantification. It performs a quasi-mapping of reads to the
    transcriptome and uses a lightweight model to estimate 
    transcript abundance. Salmon has better bias correction
    models than other tools like Kallisto and it is also
    more accurate than Kallisto.
    @Input:
        Paired-end FastQ files (scatter-per-sample)
    @Output:
        Per-sample transcript abundance estimates
    """
    input:
        r1 = join(workpath, "fastqs", "{name}.R1.fastq.gz"),
        r2 = join(workpath, "fastqs", "{name}.R2.fastq.gz"),
        index = join(workpath, "temp", "salmon_index", "seq.bin"),
    output:
        salmon = join(workpath, "counts", "transcripts", "{name}", "quant.sf"),
    params:
        rname  = "salmonquant",
        outdir = join(workpath, "counts", "transcripts", "{name}"),
        prefix = join(workpath, "temp", "salmon_index"),
        tmpdir = join(workpath, "temp"),
    resources:
        mem   = allocated("mem",  "isoformswitchanalyzer_salmon_quant", cluster),
        time  = allocated("time", "isoformswitchanalyzer_salmon_quant", cluster),
    threads: int(allocated("threads", "isoformswitchanalyzer_salmon_quant", cluster))
    container: config["images"]["isoformswitchanalyzer"]
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    export TMPDIR="${{tmp}}"

    # Quantifies transcript abundance using Salmon
    salmon quant \\
        --threads {threads} \\
        -i {params.prefix} \\
        -l A \\
        -1 {input.r1} \\
        -2 {input.r2} \\
        -o {params.outdir} \\
        --validateMappings \\
        --threads {threads} \\
        --seqBias \\
        --gcBias \\
        --posBias
    """

