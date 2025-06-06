#!/usr/bin/env Rscript

# Load necessary libraries
suppressMessages(library("leafcutter"))

# Misc helper functions
err <- function(...){cat(sprintf(...), sep='\n', file=stderr())}
fatal <- function(...) {err(...); quit(status = 1)}
timestamp <- function(..., sep = " ") {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s]", ts), ..., "\n", sep = sep)
}


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    fatal("Usage: intron_annotation.R <path_to_prepare_results_rdata> <output_intron_anno_tsv>")
}
rdata_file <- args[1]
output_tsv <- args[2]

# Load previous rdata environment
timestamp('Reading in Rdata file: ', rdata_file)
load(rdata_file)

# Print list of variables in Rdata file,
# maybe useful for debugging if something
# changes later with leafviz (i.e introns 
# variable does not exist)
timestamp('List of variables within Rdata file: ', ls())
timestamp('Here is a preview of the introns variable:')
str(introns)

# Write intron annotations to output
timestamp('Write intron annotations to output file: ', output_tsv)
write.table(introns, file=output_tsv, row.names = F, sep = "\t", quote = F)
