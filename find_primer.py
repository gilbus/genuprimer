import argparse
import configparser
import logging
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
    # setup logging and level
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    # parse all arguments, let the argparse-module do its wonderful work
    args = parse_arguments()
    logging.debug('Received arguments: {}'.format(args))
    # Parsing of the config file to get primer3-settings
    config = configparser.ConfigParser()
    config.read(args.config)
    logging.debug('Successfully read config {}'.format(args.config.name))
    # extract sequence from sequences
    sequence = parse_fasta(args.FastaFile, args.sequence)
    if not sequence:
        sys.exit("Could not find sequence with given ID-Prefix")
    # are there any settings for primer3?
    if config.has_section('primer3'):
        primer3_conf = config['primer3']
        logging.debug('Found settings for primer3, containing {}'.format(
            primer3_conf.keys()))
    else:
        primer3_conf = None
        logging.debug('Config contains no settings for primer3')

    # generate primers, dependent on sequence, primer3-configuration and names of the files containing the found primers
    generate_primer(sequence, primer3_conf, args.primerfiles)

    # is an already existing bowtie-index specified?
    if not args.index:
        # no index available, so we have to create our own one
        bowtie_index = setup_bowtie(args.FastaFile)
    else:
        bowtie_index = args.index

    # any settings for bowtie?
    if config.has_section('bowtie'):
        bowtie_conf = config['bowtie']
    else:
        bowtie_conf = None
    run_bowtie(bowtie_index, args.primerfiles, args.bowtie, bowtie_conf)


def setup_bowtie(fasta_file):
    """
    If no bowtie-index is specified we have to build it
    :type fasta_file: argparse.FileType
    """
    import re
    bowtie_index_dir = 'bowtie-index'
    # create new directory for the index, no problem if specific folder already exists
    os.makedirs(bowtie_index_dir, exist_ok=True)
    # determine name for index from name of the FASTA-file containing the sequences
    bowtie_index = "{index_dir}/{prefix}_bowtie".format(
        index_dir=bowtie_index_dir,
        prefix=
        re.split("/|\.", fasta_file.name)[-2])
    subprocess.call(["bowtie-build", fasta_file.name, bowtie_index])
    return bowtie_index


def run_bowtie(bowtie_index, files_prefix, bowtie_exec, bowtie_config):
    """
    Calls bowtie to execute the search for matches of the designed primers with other sequences.
    :param bowtie_exec: bowtie-executable, either default or defined by user
    :param bowtie_index: location of the index for bowtie
    :param files_prefix: the prefix of the files containing the primer
    :return:
    """
    # determine name of files where the previously found primers are stored
    left = "{}_left.fas".format(files_prefix)
    right = "{}_right.fas".format(files_prefix)
    # base for calling bowtie
    args = [bowtie_exec, "-k", "5000", "-S", "-f", bowtie_index, "-1",
            left, "-2", right, "--sam-nohead"]
    if bowtie_config is not None:
        if bowtie_config.getint('MaxInsSize', -1) != -1:
            args += ['-X', bowtie_config.getint('MaxInsSize')]
        if bowtie_config.getint('MinInsSize', -1) != -1:
            args += ['-I', bowtie_config.getint('MinInsSize')]
    logging.debug(args)
    res = subprocess.check_output(args)
    print(str(res).encode('UTF-8'))


def generate_primer(sequence, config, primer_file_prefix):
    """
    Calls the primer3-module with the settings and separates the results in left and right primer pairs.
    :param sequence: Sequence template for the primers
    :param config: containing the settings for primer3
    :return: The separated primer-pairs.
    """
    primer3_config_dict = {}  # type: dict
    # are any specific settings for primer3?
    if config is not None:
        # if yes, we have to parse them
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

    # set primer3-settings
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
            # only search for RIGHT or LEFT; MIDDLE would be available as well but
            # we only want the primer
            if line[1] in ['RIGHT'] and \
                    line[3] == 'SEQUENCE':
                primer_right.update({k: res[k]})
            elif line[1] in ['LEFT'] and \
                    line[3] == 'SEQUENCE':
                primer_left.update({k: res[k]})
    # write the found primer to their corresponding files
    primerfile_left = open(
        '{prefix}_left.fas'.format(prefix=primer_file_prefix), 'x')
    primerfile_right = open(
        '{prefix}_right.fas'.format(prefix=primer_file_prefix), 'x')

    # create function for sorting the prefixes
    def primer_sort(x):
        return x.split('_')[2] + x.split('_')[1]

    for k in sorted(primer_left.keys(), key=primer_sort):
        # format in FASTA-style
        line = ">{}\n{}\n\n".format(k, primer_left[k])
        primerfile_left.write(line)
    primerfile_left.close()

    for k in sorted(primer_right.keys(), key=primer_sort):
        # format in FASTA-style
        line = ">{}\n{}\n\n".format(k, primer_right[k])
        primerfile_right.write(line)
    primerfile_right.close()


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
        "--bowtie", type=str,
        default="bowtie",
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
        logging.debug('No seq_id passed, taking first sequence from {}'.format(
            fasta_file.name))
        # read one line to remove headerline from internal buffer
        fasta_file.readline()
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

    logging.debug('Successfully extracted sequence')
    return seq


if __name__ == "__main__":
    main()
