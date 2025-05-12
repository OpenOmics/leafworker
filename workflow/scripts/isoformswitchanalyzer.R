#!/usr/bin/env Rscript

########################################################
# ABOUT:
# This script identifies and quantifies isoform switches
# using the R package IsoformSwitchAnalyzeR. It requires
# the a sample sheet with sample to group information,
# a TSV file to map each sample to their salomn quant
# file, a FASTA of the transcriptome, a GTF to define
# gene models, a two groups to compare, and an output
# directory to store the results. 
########################################################


# Load the IsoformSwitchAnalyzeR library
suppressMessages(library(IsoformSwitchAnalyzeR))  # Bioconductor package 
suppressMessages(library(argparse))               # CRAN package


###################################
# Helper functions
###################################

# Create an switchAnalyzeRlist object
# by import salmon counts, design, and
# annotation information:
# https://rdrr.io/bioc/IsoformSwitchAnalyzeR/man/importRdata.html
isa_import <- function(
        sample2counts,  # TSV file to map each sample to their counts file
        output_dir,     # Output directory to write the results   
        sample_sheet,   # Sample sheet (TSV), needs sampleID and condition cols
        gtf_file,       # GTF file used for alignment and isoform quantification
        transcripts_fa, # Path to the transcriptomic fasta file
        condition_1,    # Group1 in the contrast, creates comparison g1 vs. g2
        condition_2     # Group2 in the contrast, creates comparison g1 vs. g2
    ) {

    # Import the sample sheet with
    # sample to group mappings,
    # which is used to create the
    # design matrix, must contain
    # the following columns:
    #     - 'sampleID'
    #     - 'condition'
    design_matrix <- read.table(
        sample_sheet,
        header = TRUE,
        sep = "\t",
        quote = ""
    )

    # Filter the design matrix to only
    # include the samples of interest
    design_matrix <- design_matrix[
        design_matrix$condition == condition_1 | design_matrix$condition == condition_2,
    ]

    # Create a contrast/comparison
    # data frame object, must contain
    # the following columns:
    #     - 'condition_1' ~ i.e Case
    #     - 'condition_2' ~ i.e Control (baseline group)
    comparison <- data.frame(
        condition_1 = condition_1,
        condition_2 = condition_2
    )

    # Import the per-sample salmon
    # results using sample vector,
    # this maps each sample to a
    # salmon quant.sf file
    sample_counts_location <- read.table(
        sample2counts,
        header = TRUE,
        sep = "\t",
        quote = ""
    )
    # First column is basename of sample,
    # Second columns is abs path to that
    # sample's counts file
    sample_vector        <- sample_counts_location[, 2]
    names(sample_vector) <- sample_counts_location[, 1]
    exp <- importIsoformExpression(
        sampleVector = sample_vector,
        addIsofomIdAsColumn = TRUE,
        normalizationMethod = 'TMM',
        calculateCountsFromAbundance = TRUE,
        ignoreAfterBar = TRUE,
        ignoreAfterPeriod = FALSE,
        showProgress = TRUE
    )

    # Create SwitchAnalyzeRlist object
    # from provided inputs
    isa_list <- importRdata(
        isoformCountMatrix = exp$counts,
        isoformRepExpression = exp$abundance,
        isoformExonAnnoation = gtf_file,
        designMatrix = design_matrix,
        comparisonsToMake = comparison,
        ignoreAfterBar = TRUE,
        ignoreAfterPeriod = FALSE,
        showProgress = TRUE,
        isoformNtFasta = transcripts_fa
    )


    # Filter the SwitchAnalyzeRlist,
    # extra options are current set
    # to the packages defaults,
    # using default filters
    isa_list <- preFilter(
        switchAnalyzeRlist = isa_list,
        geneExpressionCutoff = 1,
        isoformExpressionCutoff = 0,
        IFcutoff = 0.01,
        alpha = 0.05,
        dIFcutoff = 0.1,
        quiet = FALSE
    )

    return(isa_list)
}


###################################
# Parse Arguements
###################################

# Pass command line args to main
# create parser object
parser <- ArgumentParser()

# TSV file to map each sample to
# its per-sample counts files.
# Example:
# Name  File
# s1    /path/to/s1/quant.sf
# s2    /path/to/s2/quant.sf
# s3    /path/to/s3/quant.sf
# s4    /path/to/s4/quant.sf
# s5    /path/to/s5/quant.sf
# s6    /path/to/s6/quant.sf
# s7    /path/to/s7/quant.sf
# s8    /path/to/s8/quant.sf
parser$add_argument(
    "-i", "--input_sample_to_counts_file",
    help = "A TSV file containing the basename of each sample and the absolute path to its per sample counts, must contain two columns 'Name' and 'File'.",
    type = "character",
    required = TRUE
)

# Output directory to write the
# IsoformSwitchAnalyzeR results
parser$add_argument(
    "-o", "--output_directory",
    help = "Path to output directory",
    type = "character",
    required = FALSE,
    default = "IsoformSwitchAnalyzeR_results"
)

# Path to the GRCh38 plus ISO-seq
# combined GTF file
parser$add_argument(
    "-g", "--gtf_file",
    help = "Path to the GTF file used for alignment and isoform quantification",
    type = "character",
    required = TRUE
)

# Path to transcriptomic FASTA file
parser$add_argument(
    "-t", "--transcriptome_fa",
    help = "Path to the transcriptomic fasta file",
    type = "character",
    required = TRUE
)

# Sample sheet in TSV format
# containing group information
# for each samples
# Example:
# sampleID  condition   cov1
# s1        KO          B1
# s2        KO          B0
# s3        KO          B1
# s4        KO          B0
# s5        WT          B1
# s6        WT          B0
# s7        WT          B1
# s8        WT          B0
parser$add_argument(
    "-s", "--sample_sheet",
    help = "Path to the sample sheet, must contain 'sampleID' & 'condition' columns, remaining columns are treated as covariates.",
    type = "character",
    required = TRUE
)

# Group1 in a contrast to analyze,
# the resulting contrast will be
# group1 vs. group2
parser$add_argument(
    "-c1", "--condition_1",
    help = "Group1 in the contrast, the contrast will be 'group1 vs. group2'",
    type = "character",
    required = TRUE
)

# Group2 in a contrast to analyze,
# the resulting contrast will be
# group1 vs. group2
parser$add_argument(
    "-c2", "--condition_2",
    help = "Group2 in the contrast, the contrast will be 'group1 vs. group2'. This represents the baseline group in the comparison.",
    type = "character",
    required = TRUE
)

# Parse the command line arguments
args <- parser$parse_args()

# Create output directory
# if it does not exist
outdir <- args$output_directory
comparison <- paste(args$condition_1, args$condition_2, sep = "-")
prefix <- file.path(outdir, comparison)
dir.create(
    file.path(outdir),
    showWarnings = FALSE,
    recursive = TRUE
)

# Main Entry Point of Program
isa_list <- isa_import(
    sample2counts = args$input_sample_to_counts_file,
    output_dir = args$output_directory,
    sample_sheet = args$sample_sheet,
    gtf_file = args$gtf_file,
    condition_1 = args$condition_1,
    condition_2 = args$condition_2,
    transcripts_fa = args$transcriptome_fa
)

# Print summary
summary(isa_list)

# Run the main analysis, test for 
# isoform switches using DEX-seq, 
# setting reduceToSwitchingGenes
# to TRUE will cause the function
# to subset the switchAnalyzeRlist
# to the genes which each contain 
# at least one differential used 
# isoform, as indicated by the 
# alpha and dIFcutoff cutoffs
isa_list <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = isa_list,
    reduceToSwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
)

# Run alternative splicing analysis
# to quantifty different alternative
# splicing events, such as exon skipping,
# alternative 5' and 3' splice sites,
# intron retention, etc.
isa_list <- analyzeAlternativeSplicing(
    switchAnalyzeRlist = isa_list,
    onlySwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1,
    showProgress = TRUE,
    quiet = FALSE
)

# Visualize the results
# Create splicing summary plot
pdf(paste(prefix, "_splicing_summary.pdf", sep = ""))
extractSplicingSummary(
    isa_list,
    splicingToAnalyze = 'all',
    asFractionTotal = FALSE,
    onlySigIsoforms = FALSE,
    plotGenes = FALSE,
    localTheme = theme_bw(),
)
dev.off()

# Create splicing enrichment plot
pdf(paste(prefix, "_splicing_enrichment.pdf", sep = ""))
splicing_enrichment <- extractSplicingEnrichment(
    isa_list,
    splicingToAnalyze = 'all',
    alpha = 0.05,
    dIFcutoff = 0.1,
    plot = TRUE,
    localTheme = theme_bw(base_size = 14),
    minEventsForPlotting = 10,
    returnResult = TRUE,
    returnSummary = TRUE
)
dev.off()

# Create genome wide splicing plot
pdf(paste(prefix, "_genome_wide_splicing.pdf", sep = ""))
genome_wide_splicing <- extractSplicingGenomeWide(
    isa_list,
    featureToExtract = 'all',
    splicingToAnalyze = 'all',
    alpha = 0.05,
    dIFcutoff = 0.1,
    log2FCcutoff = 1,
    violinPlot = TRUE,
    alphas = c(0.05, 0.001),
    localTheme = theme_bw(),
    plot = TRUE,
    returnResult = TRUE
)
dev.off()

# Extract the top gene and isoform switches
top_gene_switches <- extractTopSwitches(
    switchAnalyzeRlist = isa_list,
    filterForConsequences = FALSE,
    extractGenes = TRUE,    # extract genes
    alpha = 0.05,
    dIFcutoff = 0.1,
    sortByQvals = TRUE,
    n = Inf
)

top_isoform_switches <- extractTopSwitches(
    switchAnalyzeRlist = isa_list,
    filterForConsequences = FALSE,
    extractGenes = FALSE,   # extract isoforms
    alpha = 0.05,
    dIFcutoff = 0.1,
    sortByQvals = TRUE,
    n = Inf
)

# Write the results to output files
# Gene switches
write.table(
    top_gene_switches,
    file = paste(prefix, "_top_gene_switches.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
# Isoform switches
write.table(
    top_isoform_switches,
    file = paste(prefix, "_top_isoform_switches.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
# Splicing enrichment
write.table(
    splicing_enrichment,
    file = paste(prefix, "_splicing_enrichment.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
# Genome wide splicing
write.table(
    genome_wide_splicing,
    file = paste(prefix, "_genome_wide_splicing.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Save all R objects to an Rdata
# file for future figures or analysis
save.image(file = paste(prefix, "_IsoformSwitchAnalyzeR.Rdata", sep = ""))

