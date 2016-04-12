import argparse
import configparser
import os
import subprocess
import sys
from argparse import RawTextHelpFormatter

import primer3


def main():
    """
    Delegates all tasks to the other functions.
    :return:
    """
    # parse all arguments, let the argparse-module do its wonderful work
    args = parse_arguments()
    # extract sequence from sequences
    sequence = parse_fasta(args.FastaFile, args.sequence)
    if not sequence:
        sys.exit("Could not find sequence with given ID-Prefix")
    primer_left, primer_right = find_primer(sequence, args.config)
    primerfile_left = open('{prefix}_left.fas'.format(prefix=args.primerfiles),
                           'x')
    primerfile_right = open(
        '{prefix}_right.fas'.format(prefix=args.primerfiles), 'x')
    for k in sorted(primer_left.keys(),
                    key=lambda x: x.split('_')[2] + x.split('_')[1]):
        line = ">{}\n{}\n\n".format(k, primer_left[k])
        primerfile_left.write(line)
    primerfile_left.close()
    for k in sorted(primer_right.keys(),
                    key=lambda x: x.split('_')[2] + x.split('_')[1]):
        line = ">{}\n{}\n\n".format(k, primer_right[k])
        primerfile_right.write(line)
    primerfile_right.close()
    bowtie_index = args.index
    if not bowtie_index:
        bowtie_index = setup_bowtie(args.FastaFile)
    run_bowtie(bowtie_index, args.primerfiles, args.bowtie)


def setup_bowtie(fasta_file):
    """
    If no bowtie-index is specified we have to build it
    :type fasta_file: argparse.FileType
    """
    import re
    bowtie_index_dir = 'bowtie-index'
    os.makedirs(bowtie_index_dir, exist_ok=True)
    bowtie_index = "{index_dir}/{prefix}_bowtie".format(
        index_dir=bowtie_index_dir,
        prefix=
        re.split("/|\.", fasta_file.name)[-2])
    subprocess.call(["bowtie-build", fasta_file.name, bowtie_index])
    return bowtie_index


def run_bowtie(bowtie_index, files_prefix, bowtie_exec):
    """
    Calls bowtie to execute the search for matches of the designed primers with other sequences.
    :param bowtie_index: location of the index for bowtie
    :param files_prefix: the prefix of the files containing the primer
    :return:
    """
    left = "{}_left.fas".format(files_prefix)
    right = "{}_right.fas".format(files_prefix)
    res = subprocess.check_output(
        [bowtie_exec.name, "-k", "5000", "-S", "-f", bowtie_index, "-1",
         left, "-2", right, "--sam-nohead"]
    )
    print(str(res))


def find_primer(sequence, configfile):
    """
    Calls the primer3-module with the settings and separates the results in left and right primer pairs.
    :param sequence: Sequence template for the primers
    :param configfile: Configfile containing the settings for primer3
    :return: The separated primer-pairs.
    """
    # Parsing of the config file to get primer3-settings
    config = configparser.ConfigParser()
    config.read(configfile.name)
    primer3_config_dict = {}  # type: dict
    for k in config['primer3'].keys():
        # Values are either ints or floats, but the primer3-module does not
        # a float where an int is required
        try:
            # if it is a int, treat it as such
            value = int(config['primer3'][k])
        except ValueError:
            # not int -> must be float
            value = float(config['primer3'][k])
        primer3_config_dict.update({str(k).upper(): value})
    # set global primer3-settings which would be reused in another run
    primer3.setP3Globals(primer3_config_dict)

    # remove any newlines or anything else like that
    sequence = sequence.replace('\n', '').replace('\r', '')
    res = primer3.bindings.designPrimers({
        'SEQUENCE_ID': 'mySequence',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': [36, 342]
    })
    primer_left = {}  # type: dict
    primer_right = {}  # type: dict
    for k in res.keys():
        line = k.split('_')
        # only applies to keys containing primer-sequences
        if len(line) >= 4:
            if line[1] in ['RIGHT'] and \
                    line[3] == 'SEQUENCE':
                primer_right.update({k: res[k]})
            elif line[1] in ['LEFT'] and \
                    line[3] == 'SEQUENCE':
                primer_left.update({k: res[k]})
    return primer_left, primer_right


def parse_arguments():
    """
    This function parses the commandline arguments via the argparse-module, which additionally generates
    a help-hook, if a parameter is passed wrong.
    :return: A parser-object containing all parsed values.
    """
    parser = argparse.ArgumentParser(
        description=
        """
        Call primer3 to generate primers and validate
        their uniqueness among other
        sequences with the help of bowtie.
        """
        , formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "FastaFile", type=argparse.FileType('r'),
        help="File containing the sequences in FASTA-Format"
    )
    parser.add_argument(
        "-s", "--sequence", type=str,
        default=None,
        help=
        """
        ID of the sequence which should be given to primer3
        (default: first sequence in FASTA-File).
        If omitted prefix-matching is used for identification
        and first hit will be used.
        """
    )
    parser.add_argument(
        "-p", "--primerfiles", type=str, default="primer3", help=
        """
        Prefix for the files where the primer will be stored. The suffixes will be 'left' and 'right'.
        Therefore the default remains, their names will 'primer_left.fas' and 'primer_right.fas'
        """
    )
    parser.add_argument(
        '-i', '--index', type=str, default=None, help=
        """
        Use existing bowtie-index. This option is directly forwarded to bowtie.
        """
    )
    parser.add_argument(
        '-c', '--config', type=argparse.FileType("r"),
        help=
        """
        Configfile with various parameters which are
        passed through to primer3. Has to include a
        'default'-section
        at top of the file and a 'primer3' with the settings.
        A default config is distributed with this programm.
        Non-working example:
        [default]
        # other settings

        [primer3]
        PRIMER_OPT_SIZE = 14,
        # more settings
        """, required=True)
    parser.add_argument(
        "--bowtie", type=argparse.FileType("r"),
        default=subprocess.check_output(["which", "bowtie"]),
        help=
        """
        The bowtie-executable if not in PATH.
        """
    )
    return parser.parse_args()


def parse_fasta(fasta_file, seq_id):
    """
    Parses the submitted FASTA-File and extracts the sequence for primer3.
    :type fasta_file: argparse.FileType
    :type seq_id: str
    :param seq_id: identifier of the sequence, for which primers should be generated
    :param fasta_file: already readable-opened file which contains all the sequences
    """
    seq = ""
    # no sequence-id specified, therefore the first one is taken
    if not seq_id:
        # fasta-style-id is found after the first sign ('>')
        seq_id = fasta_file.readline()[1:]
        for line in fasta_file:
            # read as long as no newline appears
            if line in ['\n', '\r\n']:
                break
            seq += line
    # a sequence-id has been passed to the programm
    else:
        for line in fasta_file:
            # begin of a new sequence and prefix-matching says true
            if line[0] == '>' and line.startswith(seq_id, 1):
                break
        for line in fasta_file:
            # again, read as long as no newline appears
            if line in ['\n', '\r\n']:
                break
            seq += line

    return seq


if __name__ == "__main__":
    main()
