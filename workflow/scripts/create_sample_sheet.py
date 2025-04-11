#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Author: Skyler Kuhn

# Standard Library
from __future__ import print_function
import argparse, textwrap, re, sys

# CONSTANTS
_HELP = textwrap.dedent("""
@Usage:
    $ ./create_sample_sheet.py [-h] [--version] \\
            --case-group CASE_GROUP \\
            --control-group CONTROL_GROUP \\
            --group-file GROUP_FILE

@About:
    Given a case and control group, this script will
    create a sample sheet filtered for those two
    groups. Leafcutter expects a sample sheet with
    only two groups, i.e. case vs. control, where
    the first group in the sample sheet becomes
    the baseline group. 
    The resulting sample sheet is tab-delimited,
    and it will be output to standard output.

@Required Arguments:
  --group-file   GROUPS_FILE
                   Groups file from the pipeline.
                   The second column in this file
                   contains group information.

                        @Options
  --case-group     CASE_GROUP
                    The group name for the case group.
                    This is the first group in a contrast,
                    "CASE vs. control".
  --control-group  CONTROL_GROUP
                    The group name for the control group.
                    This is the second group in a contrast,
                    "case vs. CONTROL". It represents the
                    baseline group, and it will become the
                    first group in the sample sheet.

@Optional Arguments:
  --h, --help      Shows this help message and exits.
"""
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


def clean(s, remove=['"', "'"]):
    """Cleans a string to remove any defined leading or trailing characters.
    @param s <str>:
        String to clean.
    @param remove list[<str>]:
        List of characters to remove from beginning or end of string 's'.
    @return s <str>:
        Cleaned string
    """
    for c in remove:
        s = s.strip(c)
    return s


def parse_groups(input_file, search_group, skip_header=True, delim='\t'):
    """Parses the groups file and prints lines to standard output
    with contain the group of interest.
    and their corresponding groups.
    @param input_file <file>:
        The input file to parse.
    @param search_group <str>:
        The group name to filter for.
    @param skip_header <bool>:
        Whether to skip the header (first) line or not,
        default: True
    @param delim <str>:
        Delimiter of the input_file,
        default: tab character
    @return None, prints to standard output.
    """
    with open(input_file, 'r') as input_file:
        if skip_header:
            # Ignore first line of file
            _ = next(input_file)
        for line in input_file:
            linelist = [l.strip() for l in line.split(delim)]
            try:
                group = linelist[1]    # Group column 
                if not group: continue # skip over empty string
            except IndexError:
                # No group information, skip over line
                continue
            # Check for multiple groups,
            # split on comma or semicolon
            multiple_groups = re.split(';|,', group)
            multiple_groups = [clean(g.strip()) for g in multiple_groups]
            if search_group in multiple_groups:
                # Got a hit print formatted line
                # with sample, group, and N covariates,
                # Make sure any white spaces are removed.
                covariates = "\t".join([cv.replace(" ", "") for cv in linelist[2:]])
                if covariates:
                    covariates = "\t{0}".format(covariates)
                print(
                    "{0}\t{1}{2}".format(
                        linelist[0],                    # sample column
                        search_group.replace(" ", ""),  # group column
                        covariates                      # N remaining columns (covariates)
                    )
                )


def main():
    """Collect command line args and create sample sheet."""
	# Parse command-line arguments
    # Create a top-level parser
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=_HELP,
        add_help=False
    )

    # Required Arguments
    # Groups TSV file
    parser.add_argument(
        '--group-file', required=True, help = argparse.SUPPRESS
    )

    # Case group name for filtering
    # groups tsv file
    parser.add_argument(
        '--case-group', required=True, type=str, help = argparse.SUPPRESS
    )
    # Control group name for filtering
    # groups tsv file
    parser.add_argument(
        '--control-group', required=True, type=str, help = argparse.SUPPRESS
    )

    # Optional Arguments
    # Add custom help message
    parser.add_argument(
        '-h', '--help', action='help', help=argparse.SUPPRESS
    )

    # Collect parsed arguments
    args = parser.parse_args()

    # Sanity check for usage
    if len(sys.argv) == 1:
        # Nothing was provided
        fatal('Invalid usage: create_sample_sheet.py [-h] ...')

    # Parse baseline group (i.e control group),
    # This should be the first group in the
    # leafcutter groups file.
    parse_groups(args.group_file, search_group=args.control_group)
    # Parse case group, this should be the second
    # group in the leafcutter groups file.
    parse_groups(args.group_file, search_group=args.case_group)


if __name__ == '__main__':
    # Call main method
    main()
