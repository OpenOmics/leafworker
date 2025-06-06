#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
from textwrap import dedent
import argparse, gzip, os, sys

# Constants
# Usage and help section 
_HELP = dedent("""
@Usage:
    $ ./exon_annotation.py [-h] [--version] \\
            --input GTF_FILE \\
            --output TSV_FILE 

@About:
    Given an annotation in GTF format, this script
    parses out the following information into a tab-
    separated (TSV) output file:
        • gene_id
        • gene_name
        • transcript_id
        • transcript_name
        • exon_id
        • exon_number
        • exon_start
        • exon_end
        • exon_length
    
    If the 9th column of the GTF file does not
    contain the attributes listed above, they
    will resolve to "Unknown" in the output file.

@Required:
    -i, --input GTF_FILE
        Input GTF file to parse exon information.
        This script can accept a gzip compressed
        file as input.
    -o, --output TSV_FILE
        Output TSV file with parsed exon information.

@Options:
    -h, --help     Shows this help message and exits.
    -v, --version  Prints the version and exits.

@Example:
    $ ./exon_annotation.py \\
            --input gencode.v48.annotation.gtf.gz \\
            --output exons.tsv
"""
)

# Semantic version
_VERISON = '1.0.0'


# Helper functions
def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)



def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def parse_cli_arguments():
    """Parses command line arguments and returns
    an argparse.parse_args object.
    @return <argparse.parse_args()>:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        add_help=False,
        description=_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage = argparse.SUPPRESS,
    )
    # Input GTF file to parse
    parser.add_argument(
        '-i', '--input',
        type=str, required=True,
        help=argparse.SUPPRESS
    )
    # Parsed output TSV file
    parser.add_argument(
        '-o', '--output',
        type=str, required=True,
        help=argparse.SUPPRESS
    )
    # Get version information
    parser.add_argument(
        '-v', '--version',
        action='version',
        help = argparse.SUPPRESS,
        version='%(prog)s {0}'.format(_VERISON)
    )
    # Add custom help message
    parser.add_argument(
        '-h', '--help',
        action='help',
        help=argparse.SUPPRESS
    )
    return parser.parse_args()


def stripped(s):
    """Cleans string to remove quotes
    @param s <str>:
        String to remove quotes or clean
    @return s <str>:
        Cleaned string with quotes removed
    """
    return s.strip('"').strip("'")


def parsed_attributes(attributes, parse):
    """Parses the attributes (9th column) of the GTF file
    and returns a dictionary of key-value pairs.
    @param attributes <str>:
        The attributes string from the GTF file. This should
        be the 9th column of the GTF file (1-based counting).
        It should be in the format of key "value" pairs
        separated by semicolons, where each key, value pair
        is seperated by a single space (GTF specification). 
    @param parse <list[str]>:
        List of attributes to parse from the GTF file.
        If an attribute is not found in the in the
        attributes string, it will be set to "Unknown".
    @return <dict[str]=str>:
        Dictionary containing parsed attributes
        key, value pairs.
    """
    # Setting default values for attributes
    parsed_attributes = {k: "Unknown" for k in parse}
    # Parse attributes from the GTF file
    attributes_list = attributes.split(';')
    attr_dict = {}
    for attr in attributes_list:
        if not attr.strip():
            continue # skip empty field
        # Key, value pair separated by 1 space
        key_value = attr.strip().split(' ', 1)
        if len(key_value) == 2:
            key, value = key_value
            attr_dict[key] = stripped(value)
    # Add parsed attributes to the dictionary
    for k in parse:
        if k in attr_dict:
            parsed_attributes[k] = attr_dict[k]
    return parsed_attributes


def parsed_feature(
        gtf_file,
        feature='exon',
        parse=[
            "gene_id","gene_name",
            "transcript_id","transcript_name",
            "exon_id","exon_number"
        ]
    ):
    """Parses the input GTF file and extracts lines containing
    exon information into a dictionary. Each key in the dict
    represents a column in the GTF file. 
    Here is an overview of the keys/columns:
      0. seqname: Name of the chromosome or scaffold
      1. source: Name of the program that generated this feature
      2. feature: Type of feature (e.g., exon)
      3. start: Start position of the feature in the sequence
      4. end: End position of the feature in the sequence
      5. score: Score of the feature (usually '.')
      6. strand: Strand of the feature ('+' or '-')
      7. frame: Frame of the feature (usually '.')
      8. attribute: Additional metadaata in key-value pairs
    @param gtf_file <str>:
        Path to the input GTF file to parse
    @param feature <str>:
        The type of feature to parse from the GTF file.
        Default is 'exon'.
    @param parse <list[str]>:
        List of attributes to parse from the GTF file.
        Default includes: 
            - gene_id, gene_name
            - transcript_id, transcript_name
            - exon_id, exon_number
    @return <dict[str]=str>:
        Dictionary containing parsed feature information.
        The attribute key contains a dictionary of parsed
        attributes from the 9th column of the GTF file.
    """
    # Handler for uncompressed files 
    open_func = open
    if gtf_file.endswith('.gz'):
        # Handler for gzip files
        open_func = gzip.open
    
    line_number = 0  # Used for error reporting 
    with open_func(gtf_file, 'rt') as fh:
        for line in fh:
            line_number += 1
            if line.startswith('#'):
                # Skip comment lines
                continue 
            # Split the line into columns
            tokens = line.strip().split('\t')
            if len(tokens) < 9:
                # Skip lines that do not have enough columns
                err("Warning: Skipping line {0}, insufficient # of columns (less than 9).".format(line_number, gtf_file))
                continue
            # Parse the feature type
            if tokens[2] == feature:
                yield {
                    "seqname": tokens[0],
                    "source": tokens[1],
                    "feature": tokens[2],
                    "start": int(tokens[3]),
                    "end": int(tokens[4]),
                    "score": tokens[5],
                    "strand": tokens[6],
                    "frame": tokens[7],
                    "attribute": parsed_attributes(tokens[8], parse=parse),
                    "length": int(tokens[4]) - int(tokens[3]) + 1
                }


if __name__ == '__main__':
    # Parse command line arguments
    args = parse_cli_arguments()
    
    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: {0} [-h] ...'.format(os.path.basename(sys.argv[0])))
    
    # Create output directory if
    # it does not exist
    output_dir = os.path.abspath(os.path.dirname(args.output))
    if not os.path.exists(output_dir):
        try: os.makedirs(output_dir)
        except OSError as e: 
            fatal(
                "Fatal error: Failed to create output directory: {0}\n{1}".format(
                    output_dir, e
                )
            )
    
    # Attributes (key, value pairs) to parse
    # from the 9th column of the GTF file
    PARSE_ATTRS=[
        "gene_id", "gene_name",
        "transcript_id", "transcript_name",
        "exon_id", "exon_number"
    ]
    # Location attributes of the exon
    LOC_ATTRS = ["start", "end", "length"]
    FEATURE = "exon"
    # Output TSV file handle
    with open(args.output, 'w') as out_fh:
        # Write header to the output file
        out_fh.write(
            "\t".join(PARSE_ATTRS + ["{0}_{1}".format(FEATURE, attr) for attr in LOC_ATTRS]) + "\n"
        )
        # Parse exon information from the GTF file
        for exon in parsed_feature(gtf_file=args.input, feature=FEATURE, parse=PARSE_ATTRS):
            # Get parsed attribute
            attr_dict = exon['attribute']
            # Parse key, value pairs attributes
            attr_list = [attr_dict.get(attr, 'Unknown') for attr in PARSE_ATTRS]
            # Parse exon location attributes
            loc_list = [str(exon[attr]) for attr in LOC_ATTRS]
            # Prepare output line
            output_line = attr_list + loc_list
            # Write to output file
            out_fh.write("{0}\n".format('\t'.join(output_line)))
 