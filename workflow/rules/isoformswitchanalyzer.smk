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
        --seqBias \\
        --gcBias \\
        --posBias
    """


rule isoformswitchanalyzer_salmon_matrix:
    """
    Data-processing step to get create a transcripts counts
    matrix of the estimated (raw) and normalized counts.
    @Input:
        Sample transcript abundance estimates (gather-across-all-samples)
    @Output:
        Raw transcripts counts matrix
        Normalized transcripts counts matrix
    """
    input:
        counts = expand(join(workpath, "counts", "transcripts", "{name}", "quant.sf"), name=samples),
    output:
        raw = join(workpath, "counts", "salmon.transcripts.raw_counts.tsv"),
        tpm = join(workpath, "counts", "salmon.transcripts.tpm_normalized.tsv"),
    params:
        rname  = "salmonmatrix",
        tmpdir = join(workpath, "temp"),
        script = join(workpath, "workflow", "scripts", "create_matrix.py"),
    resources:
        mem   = allocated("mem",  "isoformswitchanalyzer_salmon_matrix", cluster),
        time  = allocated("time", "isoformswitchanalyzer_salmon_matrix", cluster),
    threads: int(allocated("threads", "isoformswitchanalyzer_salmon_matrix", cluster))
    container: config["images"]["isoformswitchanalyzer"]
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in
    # mechanism for deletion on exit
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    export TMPDIR="${{tmp}}"

    # Build raw transcripts
    # counts matrix
    {params.script} \\
        --input {input.counts} \\
        --output {output.raw} \\
        --join-on Name \\
        --extract NumReads \\
        --use-parent-dir-as-sample-name \\
        --nan-values 0.000

    # Build normalized transcripts
    # counts matrix
    {params.script} \\
        --input {input.counts} \\
        --output {output.tpm} \\
        --join-on Name \\
        --extract TPM \\
        --use-parent-dir-as-sample-name \\
        --nan-values 0.000000

    # Rename Name column to transcript_id
    sed -i '1 s/^Name/transcript_id/' {output.raw}
    sed -i '1 s/^Name/transcript_id/' {output.tpm}
    """


rule isoformswitchanalyzer_mkgroups:
    """
    Data processing step to create a groups file for isoformswitchanalyzer
    differential switching script. This is a tab-delimited file containing
    at least two of the following columns:
        1.  sampleID: basename of each sample
        2.  condition: sample group inforamtion
        3+. remaining columns: covariates for correction
    @Input:
        Groups file provided to pipeline (indirect-gather-per-contrast)
    @Output:
        IsoformSwitchAnalyzeR groups file for a given comparison
        Sample vector file for a given comparison
    """
    input:
        sample_sheet = grp_file,
        # Get salmon counts for case and control groups
        case_counts = lambda w: expand(join(workpath, "counts", "transcripts", "{name}", "quant.sf"), name=group2samples[w.case]),
        cntl_counts = lambda w: expand(join(workpath, "counts", "transcripts", "{name}", "quant.sf"), name=group2samples[w.control]),
    output:
        grp = join(workpath, "differential_switching", batch_id, "{case}_vs_{control}", "groups_file.tsv"),
        vec = join(workpath, "differential_switching", batch_id, "{case}_vs_{control}", "sample_vector.tsv"),
    params:
        rname = "grpsfile",
        # Building string for cntrl
        # sample to their quant file
        ctrl2sf = lambda w: "\n".join([
            "{0}\t{1}".format(
                str(s),
                join(workpath, "counts", "transcripts", "{0}".format(s), "quant.sf"),
            ) for s in sorted(group2samples[w.control])
        ]),
        # Building string for case
        # sample to their quant file
        case2sf = lambda w: "\n".join([
            "{0}\t{1}".format(
                str(s),
                join(workpath, "counts", "transcripts", "{0}".format(s), "quant.sf")
            ) for s in sorted(group2samples[w.case])
        ]),
    resources:
        mem   = allocated("mem",  "isoformswitchanalyzer_mkgroups", cluster),
        time  = allocated("time", "isoformswitchanalyzer_mkgroups", cluster),
    threads: int(allocated("threads", "isoformswitchanalyzer_mkgroups", cluster))
    container: config["images"]["leafcutter"]
    shell: textwrap.dedent("""
    # Create groups file for a given comparison:
    # "{wildcards.case} vs. {wildcards.control}"
    # First group becomes baseline in contrast,
    # where:
    # 1st column    = sampleID
    # 2nd column    = condition
    # Nth column(s) = Covariates
    # Create header for file, IsoformSwitchAnalyzeR
    # expects the certain column names in design file
    paste \\
        <(head -1 {input.sample_sheet} | cut -f1 | sed 's/Sample/sampleID/') \\
        <(head -1 {input.sample_sheet} | cut -f2 | sed 's/Group/condition/') \\
        <(head -1 {input.sample_sheet} | cut -f3-) \\
    > {output.grp}
    # Create sample sheet for IsoformSwitchAnalyzeR
    awk -F '\\t' -v OFS='\\t' \\
        'NR!=1 && $2=="{wildcards.case}" {{print}}' \\
        {input.sample_sheet} \\
    >> {output.grp}
    awk -F '\\t' -v OFS='\\t' \\
        'NR!=1 && $2=="{wildcards.control}" {{print}}' \\
        {input.sample_sheet} \\
    >> {output.grp}

    # Create sampleVector file for mapping each
    # sample to its salmon quant.sf file, where:
    # 1st column  = Sample ID
    # 2nd column  = Absolute path to salmon quant.sf
    echo -e "Name\\tFile" > {output.vec}
    cat << EOF >> {output.vec}
    {params.case2sf}
    {params.ctrl2sf}
    EOF
    """)


rule isoformswitchanalyzer_diffswitching:
    """
    Data-processing step to find differential isoform switches using the R package,
    IsoformSwitchAnalyzer. This step uses the groups file created in the previous step and
    the sample vector file to find differential isoform switches. The groups file is a
    tab-delimited file containing at least two of the following columns:
        1.  sampleID: basename of each sample
        2.  condition: sample group inforamtion
        3+. remaining columns: covariates for correction
    The sample vector file is a tab-delimited file containing the following columns:
        1.  Name: basename of each sample
        2.  File: absolute path to the salmon quant.sf file
    @Input:
        Groups file for a given comparison (indirect-gather-per-contrast)
         Sample vector file for a given comparison
    @Output:
        Differential isoform switching results
    """
    input:
        grp = join(workpath, "differential_switching", batch_id, "{case}_vs_{control}", "groups_file.tsv"),
        vec = join(workpath, "differential_switching", batch_id, "{case}_vs_{control}", "sample_vector.tsv"),
    output:
        swt = join(workpath, "differential_switching", batch_id, "{case}_vs_{control}", "{case}-{control}_top_isoform_switches.tsv"),

    params:
        rname  = "diffswitch",
        outdir = join(workpath, "differential_switching", batch_id, "{case}_vs_{control}"),
        script = join(workpath, "workflow", "scripts", "isoformswitchanalyzer.R"),
        gtf    = gtf_file,
        transcripts = quantify_transcripts,
    resources:
        mem   = allocated("mem",  "isoformswitchanalyzer_diffswitching", cluster),
        time  = allocated("time", "isoformswitchanalyzer_diffswitching", cluster),
    threads: int(allocated("threads", "isoformswitchanalyzer_diffswitching", cluster))
    container: config["images"]["isoformswitchanalyzer"]
    shell: """
    # Unset any R-related env variables
    # that may have been inherited from
    # the user environment
    unset R_LIBS_USER; unset R_LIBS_SITE;
    # Run differential switching analysis for:
    #   {wildcards.case} vs. {wildcards.control}
    {params.script} \\
        -i {input.vec} \\
        -o {params.outdir} \\
        --sample_sheet {input.grp} \\
        --gtf_file {params.gtf} \\
        --transcriptome_fa {params.transcripts} \\
        --case_group {wildcards.case} \\
        --control_group {wildcards.control} \\
        --method saturn
    """
