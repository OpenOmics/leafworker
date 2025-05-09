#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Standard library
from __future__ import print_function
import argparse, sys, os, textwrap

# Third-party libraries
import pandas as pd

# Script metadata
__author__  = "Skyler Kuhn"
__version__ = "v0.1.0"

_help = textwrap.dedent(
"""./reorder_matrix.py:
Reorder/subset a CSV/TSV file based on column names.

@Usage:
    $ ./reorder_matrix.py [-h] [--version] \\
            --columns {{COL_X COL_Y COL_Z ...}} \\
            --input INPUT \\
            --output OUTPUT

@About:
    Subsets or reorders the columns of a CSV/TSV file
    based on a a space seperated list of their column
    names. The resulting file will have the same
    columns provided as input to the -c, --columns
    option.

@Required Arguments:
    -i, --input INPUT
                   Input CSV/TSV file to reorder or
                   subset. This file should be tab
                   or comma delimited. The delimeter
                   of the input file is determined by
                   the file extension.
    -o, --output OUTPUT
                   Output CSV/TSV file name. The output
                   file-type, and its delimiter, is
                   determined by the file's extension.
    -c, --columns {{COL_X COL_Y COL_Z ...}}
                   Space seperated list of column names
                   to reorder or subset in the input
                   file. One or more column names must
                   be provided.
@Options:
    --h, --help     Shows this help message and exits.

@Example:
  $ ./reorder_matrix.py -c Gene FC -i deg.tsv -o Gene2FC.tsv

"Everything that has a beginning has an end." 
  – The Oracle, The Matrix

@Author: {0}
@Version: {1}
""".format(__author__, __version__)
)


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


def parse_arguments():
    """Collect and parse command line arguments."""
    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: reorder_matrix.py -c <COL_X> -i <IN> -o <OUT>')

    # Parse command-line arguments
    # Create a top-level parser
    parser = argparse.ArgumentParser(
        usage = argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = _help,
        add_help=False
    )

    # Required Positional Arguments
    # List of input files
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help = argparse.SUPPRESS
    )
    # Output file name
    parser.add_argument(
        '-o',
        '--output',
        required=True,
        type=str,
        help = argparse.SUPPRESS
    )

    # Name of columns to reorder
    # or subset in the input file
    parser.add_argument(
        '-c',
        '--columns',
        required=True,
        type=str,
        nargs = '+',
        help = argparse.SUPPRESS
    )

    # Optional Arguments
    # Add custom help message
    parser.add_argument(
        '-h', '--help',
        action='help',
        help=argparse.SUPPRESS
    )

    # Collect parsed arguments
    args = parser.parse_args()

    return args


def reader(filename, subset=[], skip='#', **kwargs):
    """Reads in an MAF-like file as a dataframe. Determines the
    correct handler for reading in a given MAF file. Supports reading
    in TSV files (.tsv, .txt, .text, .vcf, or .maf), CSV files (.csv),
    and excel files (.xls, .xlsx, .xlsm, .xlsb, .odf, .ods, .odt ).
    The subset option allows a users to only select a few columns
    given a list of column names.
    @param filename <str>:
        Path of an MAF-like file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    # Get file extension
    extension = os.path.splitext(filename)[-1].lower()

    # Assign a handler to read in the file
    if extension in ['.xls', '.xlsx', '.xlsm', '.xlsb', '.odf', '.ods', '.odt']:
        # Read in as an excel file
        return excel(filename, subset, skip, **kwargs)
    elif extension in ['.csv']:
        # Read in as an CSV file
        return csv(filename, subset, skip, **kwargs)
    else:
        # Default to reading in as an TSV file
        # Tab is the normal delimeter for MAF or VCF files
        # MAF files usually have one of the following
        # extensions: '.tsv', '.txt', '.text', '.vcf', '.maf'
        return tsv(filename, subset, skip, **kwargs)


def excel(filename, subset=[], skip='#', **kwargs):
    """Reads in an excel file as a dataframe. The subset option
    allows a users to only select a few columns given a list of
    column names.
    @param filename <str>:
        Path of an EXCEL file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset:
        return pd.read_excel(filename, comment=skip, **kwargs)[subset]

    return pd.read_excel(filename, comment=skip, **kwargs)


def tsv(filename, subset=[], skip='#', **kwargs):
    """Reads in an TSV file as a dataframe. The subset option
    allows a users to only select a few columns given a list of
    column names.
    @param filename <str>:
        Path of an TSV file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset:
        return pd.read_table(filename, comment=skip, **kwargs)[subset]

    return pd.read_table(filename, comment=skip, **kwargs)


def csv(filename, subset=[], skip='#', **kwargs):
    """Reads in an CSV file as a dataframe. The subset option
    allows a users to only select a few columns given a list of
    column names.
    @param filename <str>:
         Path of an CSV file to read and parse
    @param subset list[<str>]:
        List of column names which can be used to subset the df
    @param skip <str>:
        Skips over line starting with this character
    @params kwargs <read_excel()>
        Key words to modify pandas.read_excel() function behavior
    @return <pandas dataframe>:
        dataframe with spreadsheet contents
    """
    if subset:
        return pd.read_csv(filename, comment=skip, **kwargs)[subset]

    return pd.read_csv(filename, comment=skip, **kwargs)


def writer(input_file, output_file, columns):
    """Takes a list of files and creates one excel spreadsheet.
    Each file will becomes a sheet in the spreadsheet where the
    name of the sheet is the basename of the file with the extension
    removed.
    @param input_file <str>:
        Input file to reorder/subset
    @param output_file <str>:
        Output filename
    @param columns list[<str>]:
        List of column names to reorder or subset
    """
    # Get file extension to
    # determine the correct
    # output file delimeter
    output_delimeter = '\t'
    file_extension = os.path.splitext(
        output_file
    )[-1].lower()
    if file_extension == '.csv':
        # Output file is a CSV file
        output_delimeter = ','

    # Create output directory as needed
    outdir = os.path.dirname(os.path.abspath(output_file))
    if not os.path.exists(outdir):
        # Pipeline output directory
        # does not exist on filesystem
        os.makedirs(outdir)

    # Create a spreadsheet from the contents of each file
    print('Reading in {}'.format(input_file))
    # Treat everything as a string to prevent
    # any data loss or precision
    df = reader(input_file, subset=columns, dtype=str)
    print('Writing to {}'.format(output_file))
    df.to_csv(output_file, index = False, header = True, sep = output_delimeter)


def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Input TSV file to reorder or subset
    inputs = args.input
    # Output file name
    output = args.output
    # List of columns for reordering
    # or subsetting
    columns = args.columns
    # Create XLSX file from the list
    # of input files
    writer(
        input_file=inputs,
        output_file=output,
        columns=columns,

    )


if __name__ == '__main__':
    # Call main method
    main()
