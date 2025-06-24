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
            --splicing-ann SPLICING_ANN_FILE \\
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
    -s, --splicing-ann SPLICING_ANN_FILE
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
    # Input isoform switch results file
    parser.add_argument(
        '-i', '--isoform-switch-results',
        type = lambda file: \
            check_permissions(parser, file, os.R_OK),
        required=True,
        help=argparse.SUPPRESS
    )
    # Input splicing annotation file
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


def parsed_isa_switching_genes(
        file, 
        geneid_column="gene_id",
        isoformid_column="isoform_id",
        fdr_column="isoform_switch_q_value",
        fdr_threshold=0.1
    ):
    """Parses an IsoformSwitchAnalyzeR results file
    and returns a dictionary of isoform switches
    with significant FDR values. The returned dict
    contains gene IDs as keys and a list of isoform
    IDs that have significant isoform switches as values.
    The FDR threshold is used to identify genes/isoforms
    with significant isoform switches. 
    @param file <str>:
        Path to the IsoformSwitchAnalyzeR results file.
    @param geneid_column <str>:
        Column name for gene IDs in the results file.
    @param isoformid_column <str>:
        Column name for isoform IDs in the results file.
    @param fdr_column <str>:
        Column name for FDR values in the results file.
    @param fdr_threshold <float>:
        FDR threshold to use for filtering isoform switches
    @return isoform_switches <dict[str]=list[str]>:
        Dictionary of genes with significant isoform switches,
        where keys are gene IDs and values are lists of
        isoform IDs that have significant isoform switches.
    """
    log("Starting to parse ISA results file: {0}".format(file))
    isoform_switches = {}
    # Handler for opening files, i.e.
    # uncompressed or gzip files
    open_func = gzip.open if file.endswith('.gz') else open
    line_number = 0  # Use for error reporting
    with open_func(file, 'rt') as fh:
        header = next(fh)
        col_idx = index_header(header)
        for line in fh:
            # Increment line number
            line_number += 1
            # Split the line into columns
            tokens = line.strip().split('\t')
            switch_fdr = float(tokens[col_idx[fdr_column]])
            if switch_fdr > fdr_threshold:
                # Skip lines with FDR above threshold
                continue
            # Get gene ID and isoform ID
            _k = tokens[col_idx[geneid_column]]  # _k represents gene_id
            if _k not in isoform_switches:
                isoform_switches[_k] = []
            isoform_switches[_k].append(
                tokens[col_idx[isoformid_column]]
            )
    log("Finished parsing ISA results file: {0}".format(file))
    return isoform_switches


def get_with_default(line_list, column_name_idx_dict, column_name, default_value="NA"):
    """Get a value from a list using the column name index
    dictionary. If the column name does not exist in the
    dictionary, return the default value (i.e "NA").
    @param line_list <list[str]>:
        List of values from a line in a file. This is the
        list that get are retrieving information from.
        This function is used to return a value (with a
        default value if missing) from within a list or
        dictionary comprehension.
    @param column_name_idx_dict <dict[str]=int>:
        Dictionary mapping column names to their index
    @param column_name <str>:
        Column name to look up in the dictionary
    @param default_value <str>:
        Default value to return if the column name does not
        exist in the dictionary. Defaults to "NA".
    @return value <str>:
        Value from the list at the index of the column name,
        or the default value if the column name does not exist.
    """
    # Default value to return if column name DNE
    parsed_value = default_value
    # Try to get the index of the column name
    # This can result in a KeyError/IndexError
    # if the column name does not exist in the
    # dict.
    if column_name in column_name_idx_dict:
        # Get the index of the column name
        # from the dictionary.
        # This will raise a KeyError if the
        # column name does not exist in the dict.
        list_idx = column_name_idx_dict[column_name]
        if list_idx < len(line_list):
            # If the index is within the bounds of the list,
            # return the value at that index.
            parsed_value = line_list[list_idx]
            # Remove quotes from the value
            parsed_value = stripped(parsed_value)
    return parsed_value


def index_file(
        file, key = "gene_id",
        single_value_columns=["gene_name"],
        multi_value_columns=["transcript_id", "transcript_name"]
    ):
    """Parses and indexes a file into a dictionary for quick
    lookups later. This file is indexed by keys representing
    columns defined via single_value_columns and multi_value_columns.
    If a columns is listed in the multi_value_columns list, it will
    return a list of values if there is a 1:M relationship.
    @param file <str>:
        File to parse and index. Must contain a header with
        the columns listed in keys and values. The index of
        these columns will be automatically resolved.
    @param key <str>:
        Key to use for the first level of the index. Defaults
        to "gene_id". This is the first key in the nested 
        dictionary.
    @param single_value_columns <list[str]>:
        List of columns that contain single values per key.
        These are the second keys in the nested dictionary.
        This returns a single value as a string.
    @param multi_value_columns <list[str]>:
        List of columns that contain multiple values per key.
        These values are stored as a list in the nested dictionary
        to accomodate 1:M relationship between the key and the
        values in this list.
    @return file_idx <dict[str]=str>:
        Nested dictionary where,
            • 1st_key = gene_id
            • 2nd_key = one of the following:
                gene_name|transcript_id|transcript_name
    @note: file_idx[transcript_id][2nd_ley] returns a
    either a string or a list, depending on whether the
    column is single_value_columns or multi_value_columns
    """
    log("Started indexing input file: {0}".format(file))
    file_idx = {}
    # Handler for opening files, i.e.
    # uncompressed or gzip files
    open_func = gzip.open if file.endswith('.gz') else open
    line_number = 0  # Used for error reporting
    with open_func(file, 'rt') as fh:
        header = next(fh)
        col_idx = index_header(header)
        for line in fh:
            # Increment line number
            line_number += 1
            # Split the line into columns
            tokens = line.strip().split('\t')
            _k = tokens[col_idx[key]]  # _k represents gene_id
            _v = {v: get_with_default(tokens,col_idx,v) for v in single_value_columns}
            if _k not in file_idx:
                # Create a new entry in the index
                # with the single value columns
                # this prevents it from being
                # overwritten in every iteration.
                # The values in _v should be the
                # same on every line for each
                # transcript_id.
                file_idx[_k] = _v
            # Parse multi-value columns
            for multi_v_key in multi_value_columns:
                if multi_v_key not in file_idx[_k]:
                    file_idx[_k][multi_v_key] = [] # init as empty list
                # Append the multi-value columns to the index
                file_idx[_k][multi_v_key].append(get_with_default(tokens,col_idx,multi_v_key))
    log("Finished indexing input file: {0} ({1} lines)".format(file, line_number))
    return file_idx


def fasta(filename):
    """
    Reads in a FASTA file and yields each of its entries.
    The generator yields each sequence identifier and its 
    corresponding sequence to ensure a low memory profile. 
    If a sequence occurs over multiple lines, the yielded 
    sequence is concatenated.
     @param filename <str>:
        Path of FASTA file to read and parse
    @yield chrom, sequence <str>, <str>:
        Yields each seq id and seq in the FASTA file
    """
    log("Started parsing FASTA file: {0}".format(filename))
    open_func = gzip.open if filename.endswith('.gz') else open
    with open_func(filename, 'rt') as file:
        sequence, chrom = '', ''
        for line in file:
            line = line.strip()
            if line.startswith('>') and sequence:
                # base case for additional entries
                yield chrom, sequence
                chrom = line[1:] # remove the > symbol
                sequence = ''
            elif line.startswith('>'):
                # base case for first entry in fasta file
                chrom = line[1:] # remove the > symbol
            else:
                # concatenate multi-line sequences
                sequence += line
        else:
            yield chrom, sequence
    log("Finished parsing FASTA file: {0}".format(filename))


if __name__ == '__main__':
    # Parse command line arguments
    args = parse_cli_arguments()

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

    # Parse the isoform switch results file
    # where keys are gene IDs and values
    # are lists of isoform IDs that have
    # significant isoform switches
    top_isa_switches = parsed_isa_switching_genes(
        args.isoform_switch_results,
        geneid_column="gene_id",
        isoformid_column="isoform_id",
        fdr_column="isoform_switch_q_value",
        fdr_threshold=args.fdr_filter
    )
    log("Found {0} genes with significant isoform switches.".format(len(top_isa_switches)))
    # Parse the splicing annotation file
    gene2isoform_dict = index_file(
        args.splicing_ann,
        key="gene_id",
        single_value_columns=["gene_name"],
        multi_value_columns=["transcript_id", "transcript_name"]
    )
    log('Indexed transcript information for {0} genes'.format(len(gene2isoform_dict)))
    # Create an index of transript_id
    # to gene_id mappings
    transcript2gene_dict = index_file(
        args.splicing_ann,
        key="transcript_id",
        single_value_columns=["gene_id"],
        multi_value_columns=[]
    )
    log('Indexed gene information for {0} transcripts'.format(len(transcript2gene_dict))) 
    # Parse the transcripts FASTA file,
    # where keys are transcript IDs and values
    # are the corresponding concatnated sequences
    transcript_sequences = {}
    for chrom, sequence in fasta(args.transcripts_fasta):
        # Extract transcript ID from the FASTA header
        # using the delimiter specified by the user
        transcript_id = chrom.split(args.fasta_delim)[0]
        transcript_sequences[transcript_id] = sequence
        transcript_sequences[transcript_id] = sequence
    log("Parsed {0} transcripts from FASTA file.".format(len(transcript_sequences)))
    # Create a dictionary to group the
    # isoform switches by gene_id, ISA
    # has a tendency to report the gene_id
    # as the gene_name, meaning that gene_id
    # and gene_name are often the same.
    log("Started grouping isoform switches by gene_id.")
    top_isa_switches_grouped_by_gene = {}
    for gene_id, transcripts_list in top_isa_switches.items():
        # Group the isoform switches by gene_id
        for transcript_id in transcripts_list:
            _converted_gene_id = transcript2gene_dict.get(transcript_id, {}).get("gene_id", "NA")
            if _converted_gene_id == "NA":
                log("Warning: Transcript ID '{0}' does not have a corresponding gene_id, skipping...".format(transcript_id))
                continue
            if _converted_gene_id not in top_isa_switches_grouped_by_gene:
                top_isa_switches_grouped_by_gene[_converted_gene_id] = []
            top_isa_switches_grouped_by_gene[_converted_gene_id].append(transcript_id)
    log("Finished grouping isoform switches by {0} gene_id.".format(len(top_isa_switches_grouped_by_gene)))
    # Write the isoform sequences to the output file
    ofh_seq_number = 0
    with open(args.output, "w") as ofh:
        for gene_id in top_isa_switches_grouped_by_gene:
            # Get the the transcript_ids for
            # all the transcripts belonging to
            # the gene_id that has a significant
            # isoform switch
            isoform_list = gene2isoform_dict.get(gene_id, {}).get("transcript_id", [])
            gene_name    = gene2isoform_dict.get(gene_id, {}).get("gene_name", "NA")
            # Write each isoform sequence to the output file
            for isoform_id in isoform_list:
                if isoform_id not in transcript_sequences:
                    log("Warning: Transcript ID '{0}' not found in FASTA file, skipping.".format(isoform_id))
                    continue
                sequence = transcript_sequences[isoform_id]
                # Output FASTA file format:
                # >transcript_id|gene_id|gene_name
                # ATGTTTACAGTACCCCAA
                ofh_seq_number += 1
                ofh.write(
                    ">{0}{1}{2}{1}{3}\n{4}\n".format(
                        isoform_id, args.fasta_delim, gene_id, gene_name, sequence
                    )
                )
    log("Finished writing {0} isoform sequences across to output file: {1}".format(ofh_seq_number, args.output))
    log("Finished running isoform sequences script!")
