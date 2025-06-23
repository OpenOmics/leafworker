#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
from datetime import datetime
from textwrap import dedent
import argparse, gzip, os, sys

# Constants
# Usage and help section
_HELP = dedent("""
@Usage:
    $ ./isoform_sequences.py [-h] [--version] \\
            [--fdr-filter FDR_FILTER] \\
            [--fasta-delim FASTA_DELIM] \\
            --isoform-switch-results ISOFORM_SWITCH_RESULTS \\
            --splicing-ann SPLICIN_ANN_FILE \\
            --transcripts-fasta TRANSCRIPTS_FASTA \\
            --output OUTPUT_FILE
@About:
    Given an IsoformSwitchAnalyzeR results file,
    a splicing annptation file, and a FASTA file
    containing transcript sequences, this script
    will parse the splicing annotation file and
    the FASTA file to generate a new FASTA file
    containing the sequences of each transcript
    with a significant isoform switch.
    This script will not only output the sequence
    of the transcript containing the isoform switch,
    but also the sequence of all the isoforms that
    belong to the same gene as the transcript
    containing the isoform switch. This allows
    for the user to see the context of the isoform
    switch and how it relates to the other isoforms
    of the same gene.

@Required:
    -i, -isoform-switch-results ISOFORM_SWITCH_RESULTS
        Input isoform switch results file. This
        file should be a tab-separated file and
        contain the following columns:
            • isoform_id
            • gene_id
            • isoform_switch_q_value
    -s, --splicing-ann SPLICIN_ANN_FILE
        Splicing annotation file. This file should
        be a tab-separated file containing the
        following columns:
            • transcript_id
            • gene_id
    -t, --transcripts-fasta TRANSCRIPTS_FASTA
        FASTA file containing transcript sequences.
        This file should contain the sequences
        of all transcripts in the splicing
        annotation file. The sequences should
        be in the FASTA format, where each
        transcript is represented by a header
        line starting with ">" followed by the
        transcript ID, and the sequence on the
        next line. The header line should
        contain the transcript ID as the first
        field, followed by any additional
        information separated by the FASTA_DELIM
        character. The FASTA_DELIM character
        defaults to "|", but can be changed
        using the --fasta-delim option.
    -o, --output OUTPUT_FILE
        Output file with transcript containing
        isoform sequences. This file will be
        a FASTA file containing the sequences
        of the transcripts that contain an
        isoform switch. The sequences will be
        the sequences of the transcripts that
        contain the isoform switch, as well as
        the sequences of all the isoforms that
        belong to the same gene as the transcript
        containing the significant isoform switch.
@Options:
    -f, --fdr-filter FDR_THRESHOLD
        FDR threshold to use for filtering isoform
        switches. This is the q-value threshold
        used to determine if an isoform switch is
        significant. Default: "0.1".
    -d, --fasta-delim FASTA_DELIM
        Delimiter character used in the FASTA
        file to separate the transcript ID from
        any additional information in the header.
        Default is "|". This character should not
        be present in the transcript ID. This char
        is used to parse the transcript ID from
        the FASTA header line. The transcript ID
        is the first field in the header line,
        followed by any additional information
        separated by the FASTA_DELIM character.
        Default: "|".
    -h, --help
        Shows help message and exits.
    -v, --version
        Prints the version and exits.

@Example:
    $ ./isoform_sequences.py \\
        -i isa_top_isoform_switches.tsv \\
        -s splicing_annotation.tsv \\
        -t transcripts.fasta \\
        -o isa_top_isoform_sequences.fa \\
        --fdr-filter 0.1 \\
        --fasta-delim "|"
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


def timestamp(format="%Y-%m-%d %H:%M:%S"):
    """Returns a formatted timestamp string
    for the current time.
    @param format <str>:
        Format string for the timestamp, default:
        "%Y-%m-%d %H:%M:%S" which is equivalent to
        "2023-10-01 12:00:00" for example.
    @return <str>:
        Formatted timestamp string, i.e. "2023-10-01 12:00:00"
    """
    return datetime.now().strftime(format)


def log(*message):
    """Logs a message to standard output with a timestamp.
    @param message <any>:
        Values printed to log
    """
    print("[{0}] {1}".format(
        timestamp(),
        " ".join([str(m) for m in message]))
    )


def check_permissions(parser, path, *args, **kwargs):
    """Checks permissions using os.access() to see the
    user is authorized to access a file/directory. Checks
    for existence, read, write and execute via args:
        • os.F_OK (tests existence)
        • os.R_OK (tests read)
        • os.W_OK (tests write)
        • os.X_OK (tests exec)
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @param args <any>:
        Positional args to pass to os.access()
    @param kwargs <any>:
        Named kwargs to pass to os.access()
    @return path <str>:
        Returns absolute path if it exists and the
        checked permssions are setup are correct.
    """
    if not os.path.exists(path):
        parser.error(
            "Path '{}' does not exists! Failed to provide vaild input.".format(path)
        )
    if not os.access(path, *args, **kwargs):
        parser.error(
            "Path '{}' exists, but cannot read path due to permissions!".format(path)
        )
    return os.path.abspath(path)


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
    # Isoform switch results file
    parser.add_argument(
        '-i', '--isoform-switch-results',
        type = lambda file: \
            check_permissions(parser, file, os.R_OK),
        required=True,
        help=argparse.SUPPRESS
    )
    # Splicing annotation file
    parser.add_argument(
        '-s', '--splicing-ann',
        type = lambda file: \
            check_permissions(parser, file, os.R_OK),
        required=True,
        help=argparse.SUPPRESS
    )
    # Input transcripts sequences FASTA file
    parser.add_argument(
        '-t', '--transcripts-fasta',
        type = lambda file: \
            check_permissions(parser, file, os.R_OK),
        required=True,
        help=argparse.SUPPRESS
    )
    # Output FASTA file of isoform sequences
    parser.add_argument(
        '-o', '--output',
        type = str,
        required=True,
        help=argparse.SUPPRESS
    )
    # FDR filtering threshold
    parser.add_argument(
        '-f', '--fdr-filter',
        type=float,
        required=False,
        help=argparse.SUPPRESS,
        default=0.1
    )
    # Input transcripts FASTA delimiter,
    # used to parse the transcript ID from
    # the FASTA header line of the file
    # provided to --transcripts-fasta
    parser.add_argument(
        '-d', '--fasta-delim',
        type=str,
        required=False,
        help=argparse.SUPPRESS,
        default="|"
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


def index_header(file_header):
    """Returns the index of each column_name
    as a dictionary.
    @param file_header <str>:
        First line of a file, containing column names
    @return idx <dict[str]=int>:
        Column name to index dictionary
    """
    idx = {}
    tokens = [
        stripped(c.strip()) \
            for c in file_header.strip().split('\t')
    ]
    # Create column name to index mapping
    for i,c in enumerate(tokens):
        idx[c]=i
    return idx


if __name__ == '__main__':
    # Parse command line arguments
    args = parse_cli_arguments()

    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: {0} [-h] ...'.format(os.path.basename(sys.argv[0])))

    log("Running isoform sequences script with the following options: ", args)
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

    log("Finished running isoform sequences script!")