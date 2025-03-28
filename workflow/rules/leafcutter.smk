# Functions and rules for running leafcutter
# Local imports
from scripts.common import (
    allocated
)


# Data processing leafcutter rules
rule leafcutter_bam2bed:
    """
    Data-processing step to keep any aligned spliced RNA reads from the input BAM
    file based on certain criteria: 
      1. Checks if the read contains a splice junction (indicated by "N"
         in the CIGAR string).
      2. Ensures the read has a quality score of at least 10.
      3. Validates that the length of the sequence mapped before the splice 
         junction (intron) is greater than 50 and that at least 6 nucleotides
         map into each exon.
    @Input:
        Input BAM file (scatter)
    @Output:
        Temp BED file of alignments meeting the criteria above
    """
    input:
        bam   = join(workpath, "{name}.bam"),
    output:
        bed   = join(workpath, "temp", "{name}.bed"),
    params:
        rname = "bam2bed",
    resources:
        mem   = allocated("mem",  "leafcutter_bam2bed", cluster),
        time  = allocated("time", "leafcutter_bam2bed", cluster),
    threads: int(allocated("threads", "leafcutter_bam2bed", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Filter alignments and create
    # intermediate bed file
    samtools view {input.bam} \\
        | filter_cs.py \\
        | sam2bed.pl \\
            --use-RNA-strand \\
            /dev/stdin \\
            {output.bed}
    """


rule leafcutter_bed2junc:
    """
    Data-processing to convert BED file in a junction file. It is worth
    noting the docs point to another tool for creating this file, regtools:
    $ regtools junctions extract -a 8 -m 50 -M 500000 $bam -o $junc
    However, the python_wrapper.sh script provided with the example_data
    use the method within this snakefile. It maybe worth switching to the
    regtools method in the future.
    @Input:
        BED file (scatter)
    @Output:
        Junction file,
        Temp file containing pointer to junction file
    """
    input:
        bed   = join(workpath, "temp", "{name}.bed"),
    output:
        jnt   = join(workpath, "junctions", "{name}.junc"),
        tmp   = join(workpath, "temp", "{name}.tmp"),
    params:
        rname = "bed2junc",
    resources:
        mem   = allocated("mem",  "leafcutter_bed2junc", cluster),
        time  = allocated("time", "leafcutter_bed2junc", cluster),
    threads: int(allocated("threads", "leafcutter_bed2junc", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Create a junction file 
    # from the BED file 
    bed2junc.pl \\
        {input.bed} \\
        {output.jnt}
    # Create tmp file with
    # pointer to junction 
    # file used for gather
    # step later
    echo '{output.jnt}' > {output.tmp}
    """


rule leafcutter_gatherjuncs:
    """
    Gather step to create a single file containing the paths to each of
    the per-sample junction files. This file is needed for leafcutter's 
    clustering step.
    @Input:
        Temp file containing pointer to junction file (gather-across-all-samples)
    @Output:
        All samples junction file
    """
    input:
        tmps  = expand(join(workpath, "temp", "{name}.tmp"), name=samples),
    output:
        jnt   = join(workpath, "temp", "junction_files.txt"),
    params:
        rname = "gthrjuncs",
    resources:
        mem   = allocated("mem",  "leafcutter_gatherjuncs", cluster),
        time  = allocated("time", "leafcutter_gatherjuncs", cluster),
    threads: int(allocated("threads", "leafcutter_gatherjuncs", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Create single file with pointing
    # to all junction files
    cat {input.tmps} \\
    > {output.jnt} 
    """


rule leafcutter_clusterjuncs:
    """
    Data processing step to pool, sort, refine, and cluster intron junctions.
    @Input:
        All samples junction file (indirect-gather-singleton)
    @Output:
        Junctions file cluster across all samples
    """
    input:
        jnt   = join(workpath, "temp", "junction_files.txt"),
    output:
        cnt   = join(workpath, "junctions", "leafcutter_perind.counts.gz"),
    params:
        rname = "clstjuncs",
        outd  = join(workpath, "junctions"),
    resources:
        mem   = allocated("mem",  "leafcutter_clusterjuncs", cluster),
        time  = allocated("time", "leafcutter_clusterjuncs", cluster),
    threads: int(allocated("threads", "leafcutter_clusterjuncs", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Cluster intron junctions
    # across all samples
    leafcutter_cluster.py \\
        -j {input.jnt} \\
        -m 30 \\
        -o leafcutter \\
        -l 100000 \\
        -r {params.outd} \\
        -p 0.001
    """
