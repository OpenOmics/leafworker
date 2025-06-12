# Functions and rules for running leafcutter
# Standard Library
import textwrap
# Local imports
from scripts.common import (
    allocated
)


# Data processing leafcutter rules
rule leafcutter_gtf2exons:
    """
    Data-processing step to create a reference file for leafcutter's 
    differential splicing script. This file is a tab-delimited file
    containing information about exons. The exons file should have
    the following columns: chr, start, end, strand, gene_name. It is
    used by leafcutter to annotate clusters to genes.
    @Input:
        Input GTF file (singleton)
    @Output:
        Exons TSV file
    """
    input:
        gtf   = gtf_file
    output:
        exons = join(workpath, "temp", "exons.tsv.gz"),
        exon_ann   = join(workpath, "temp", "exon_annotation.tsv"),
        splice_ann = join(workpath, "temp", "splicing_annotation.tsv"),
    params:
        rname = "gtf2exons",
        ex_script = join(workpath, "workflow", "scripts", "exon_annotation.py"),
        sp_script = join(workpath, "workflow", "scripts", "splicing_annotation.py"),
    resources:
        mem   = allocated("mem",  "leafcutter_gtf2exons", cluster),
        time  = allocated("time", "leafcutter_gtf2exons", cluster),
    threads: int(allocated("threads", "leafcutter_gtf2exons", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Unset any R-related env variables
    # that may have been inherited from
    # the user environment
    unset R_LIBS_USER; unset R_LIBS_SITE;
    # Get information for each exon, i.e
    # chr, start, end, strand, gene_name
    gtf_to_exons.R \\
        {input.gtf} \\
        {output.exons}
    # Parse and write exon information to
    # an output file for adding to the 
    # effect sizes output file. This file
    # contains the following information:
    # gene_id, gene_name, transcript_id, 
    # transcript_name, exon_id, exon_number,
    # exon_start, exon_end, exon_length
    echo 'Getting detailed exon information...'
    {params.ex_script} \\
        --input {input.gtf} \\
        --output {output.exon_ann}
    # Get splicing annotation for each
    # transcript, contains information
    # about the location, number, and
    # and exon_id of each exon that
    # makes up the transcript
    echo 'Getting splicing annotation...'
    {params.sp_script} \\
        --exon-ann {output.exon_ann} \\
        --output {output.splice_ann}
    """


rule leafcutter_gtf2leafviz:
    """
    Data-processing step to create a reference file for leafcutter's 
    script to visualize the differential splicing results. The shiny
    app includes functionality to label detected introns as annotated
    or cryptic. This step processes a given GTF to generate lists of
    exons, introns and splice sites.
    @Input:
        Input GTF file (singleton)
    @Output:
        All introns TSV file,
        All exons TSV file,
        5-prime introns TSV file,
        3-prime introns TSV file
    """
    input:
        gtf   = gtf_file
    output:
        introns = join(workpath, "temp", "annotation_all_introns.bed.gz"),
        exons   = join(workpath, "temp", "annotation_all_exons.txt.gz"),
        start   = join(workpath, "temp", "annotation_fiveprime.bed.gz"),
        end     = join(workpath, "temp", "annotation_threeprime.bed.gz"),
    params:
        rname  = "gtf2leafviz",
        prefix = join(workpath, "temp", "annotation"),
        tmpdir = join(workpath, "temp"),
    resources:
        mem   = allocated("mem",  "leafcutter_gtf2leafviz", cluster),
        time  = allocated("time", "leafcutter_gtf2leafviz", cluster),
    threads: int(allocated("threads", "leafcutter_gtf2leafviz", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    export TMPDIR="${{tmp}}"

    # Get detailed information about
    # exons, introns, and their start
    # and end junction positions
    gtf2leafcutter.pl \\
        -o {params.prefix} \\
        {input.gtf}
    """


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
        num   = join(workpath, "junctions", "leafcutter_perind_numers.counts.gz"),
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


rule leafcutter_clustergenes:
    """
    Data processing step to associate splicing clusters to genes.
    @Input:
        Junctions file cluster across all samples (indirect-gather-singleton)
    @Output:
        Text file containing gene splicing clusters 
    """
    input:
        cnt   = join(workpath, "junctions", "leafcutter_perind.counts.gz"),
        gtf   = gtf_file
    output:
        txt   = join(workpath, "junctions", "leafcutter.clu2gene.txt"),
    params:
        rname = "clstgene",
        outd  = join(workpath, "junctions"),
    resources:
        mem   = allocated("mem",  "leafcutter_clustergenes", cluster),
        time  = allocated("time", "leafcutter_clustergenes", cluster),
    threads: int(allocated("threads", "leafcutter_clustergenes", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Get gene splicing clusters
    get_cluster_gene.py \\
        {input.gtf} \\
        {input.cnt}
    """


rule leafcutter_mkgroups:
    """
    Data processing step to create a groups file for leafcutter's 
    differential splicing script. This is a tab-delimited file
    containing two columns: the path to a sample's bam file and
    it's group assignment. The first group in the file  will
    represent the baseline group in a contrast. 
    @Input:
        Junctions file cluster across all samples (indirect-gather-per-contrast)
    @Output:
        Groups file for a given comparison
    """
    input:
        num = join(workpath, "junctions", "leafcutter_perind_numers.counts.gz"),
    output:
        grp = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "groups_file.tsv"),
    params:
        rname = "grpsfile",
        # Building string for cntrl
        # sample to their group
        ctrl2grp = lambda w: "\n".join([
            "{0}\t{1}".format(
                str(s), str(w.control)
            ) for s in group2samples[w.control]
        ]),
        # Building string for case
        # sample to their group
        case2grp = lambda w: "\n".join([
            "{0}\t{1}".format(
                str(s), str(w.case)
            ) for s in group2samples[w.case]
        ]),
        mkgrp_script = join(workpath, "workflow", "scripts", "create_sample_sheet.py"),
        groups_file  = grp_file,
    resources:
        mem   = allocated("mem",  "leafcutter_groupsfile", cluster),
        time  = allocated("time", "leafcutter_groupsfile", cluster),
    threads: int(allocated("threads", "leafcutter_groupsfile", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Create groups file for a given comparison:
    # "{wildcards.case} vs. {wildcards.control}"
    # First group becomes baseline in contrast,
    # where:
    # 1st column    = Sample
    # 2nd column    = Group
    # Nth column(s) = Covariates
    {params.mkgrp_script} \\
        --case-group {wildcards.case} \\
        --control-group {wildcards.control} \\
        --group-file {params.groups_file} \\
    > {output.grp}
    """


rule leafcutter_diffsplicing:
    """
    Data-processing step to perform differential splicing analysis using 
    leafcutter's differential splicing script.
    @Input:
        Groups file for a given comparison (indirect-gather-per-contrast)
    @Output:
        Differential splicing results
    """
    input:
        exn = join(workpath, "temp", "exons.tsv.gz"),
        num = join(workpath, "junctions", "leafcutter_perind_numers.counts.gz"),
        grp = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "groups_file.tsv"),
    output:
        res = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "diff_splicing_cluster_significance.txt"),
        eff = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "diff_splicing_effect_sizes.txt"),
    params:
        rname  = "diffsplice",
        prefix = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "diff_splicing"),
        # Resolve min samples per intron
        # to 3 for the analysis unless 
        # there is a group with 2 samples
        # in which case set to 2
        min_samples = lambda w: 2 if min(
            len(group2samples[w.case]), len(group2samples[w.control])
        ) < 3 else 3, 
    resources:
        mem   = allocated("mem",  "leafcutter_diffsplicing", cluster),
        time  = allocated("time", "leafcutter_diffsplicing", cluster),
    threads: int(allocated("threads", "leafcutter_diffsplicing", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Unset any R-related env variables
    # that may have been inherited from
    # the user environment
    unset R_LIBS_USER; unset R_LIBS_SITE;
    # Run differential splicing analysis for:
    #   {wildcards.case} vs. {wildcards.control}
    leafcutter_ds.R \\
        --exon_file {input.exn} \\
        --num_threads {threads} \\
        --min_samples_per_intron {params.min_samples} \\
        --timeout 1200 \\
        --seed 12345 \\
        --output_prefix {params.prefix} \\
        {input.num} \\
        {input.grp}
    """


rule leafcutter_prepleafviz:
    """
    Data-processing step to create the Rdata file for running leafcutter's
    shiny application, i.e leafviz. This shiny application can be used to
    interactively explore the differential splicing results. This step will
    generate a single .RData file with everything the shiny app needs to run.
    @Input:
        Counts file provided to leafcutter_ds.R,
        Groups file provided to leafcutter_ds.R,
        Cluster significant table output by leafcutter_ds.R,
        Per-junction effect sizes table output by leafcutter_ds.R,
        All introns TSV file,

    @Output:
        Rdata file to input to leafviz shiny app
    """
    input:
        num = join(workpath, "junctions", "leafcutter_perind_numers.counts.gz"),
        grp = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "groups_file.tsv"),
        res = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "diff_splicing_cluster_significance.txt"),
        eff = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "diff_splicing_effect_sizes.txt"),
        itn = join(workpath, "temp", "annotation_all_introns.bed.gz"),
    output:
        rda = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "{case}_vs_{control}_leafviz.Rdata"),
        itn = join(workpath, "differential_splicing", batch_id, "{case}_vs_{control}", "{case}_vs_{control}_intron_annotation.tsv"),
    params:
        rname  = "prleafviz",
        prefix = join(workpath, "temp", "annotation"),
        rscript = join(workpath, "workflow", "scripts", "intron_annotation.R"),
    resources:
        mem   = allocated("mem",  "leafcutter_prepleafviz", cluster),
        time  = allocated("time", "leafcutter_prepleafviz", cluster),
    threads: int(allocated("threads", "leafcutter_prepleafviz", cluster))
    container: config["images"]["leafcutter"]
    shell: """
    # Unset any R-related env variables
    # that may have been inherited from
    # the user environment
    unset R_LIBS_USER; unset R_LIBS_SITE;
    # Create input leafviz Rdata file for:
    # "{wildcards.case} vs. {wildcards.control}"
    prepare_results.R \\
        --FDR 0.1 \\
        --meta_data_file {input.grp} \\
        --code "{wildcards.case}_vs_{wildcards.control}" \\
        {input.num} \\
        {input.res} \\
        {input.eff} \\
        {params.prefix} \\
        -o {output.rda}
    # Write the "intron" annotation table 
    # within the Rdata object to an output
    # file for adding more information to
    # the effect sizes output file
    {params.rscript} \\
        "{output.rda}" \\
        "{output.itn}"
    """
