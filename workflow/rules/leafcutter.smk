# Functions and rules for running leafcutter
# Local imports
from scripts.common import (
    allocated
)


# Data processing leafcutter rules
rule bam2bed:
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
        BED file of alignments meeting the criteria above
    """
    input:
        bam=join(workpath,"{name}.bam"),
    output:
        bed=join(workpath,"temp","{name}.bed"),
    params:
        rname="bam2bed",
    resources:
        mem       = allocated("mem",  "bam2bed", cluster),
        time      = allocated("time", "bam2bed", cluster),
    threads: int(allocated("threads", "bam2bed", cluster))
    container:  config["images"]["leafcutter"]
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


rule bed2junc:
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
        Junction file
    """
    input:
        bed=join(workpath,"temp","{name}.bed"),
    output:
        jnt=join(workpath,"temp","{name}.junc"),
    params:
        rname="bed2junc",
    resources:
        mem    = allocated("mem",  "bed2junc", cluster),
        time   = allocated("time", "bed2junc", cluster),
    threads: int(allocated("threads", "bed2junc", cluster))
    container:  config["images"]["leafcutter"]
    shell: """
    # Create a junction file 
    # from the BED file 
    bed2junc.pl \\
        {input.bed} \\
        {output.jnt}
    """
