import argparse
import ast
import configparser
import logging
import os
import subprocess
import sys
from argparse import RawTextHelpFormatter

import primer3

# mapping table for keys used in config-file and actual bowtie settings
bowtie_options_translation = {
    'MaxInsSize'.lower(): '-X',
    'MinInsSize'.lower(): '-I',
    'MaxMismatch'.lower(): '-v'
}

# 'lower()' is called since the values, once parsed from the config, are stored
# in lowercase. But it is nicer to have the same values written in the same way

# checking whether submitted values are valid and explaining messages in case
# they are not
bowtie_valid_setting = {
    'MaxMismatch'.lower(): (
        lambda x: 0 <= x <= 3, 'Value must hold 0 <= v <= 3')
}


def main():
    """
    Delegates all tasks to the other functions.
    :return:
    """
    # parse all arguments, let the argparse-module do its wonderful work
    args = parse_arguments()
    # setup logging and level as well as colors for logging
    # colors are only working on linux
    logging.addLevelName(logging.WARNING,
                         # format red
                         "\033[1;31m%s\033[1;0m" % logging.getLevelName(
                             logging.WARNING))
    logging.addLevelName(logging.INFO,
                         # format blue
                         "\033[1;34m%s\033[1;0m" % logging.getLevelName(
                             logging.INFO))
    logging.addLevelName(logging.DEBUG,
                         # format green
                         "\033[1;32m%s\033[1;0m" % logging.getLevelName(
                             logging.DEBUG))
    logging.basicConfig(format='%(levelname)s:%(message)s', level=args.loglevel)
    logging.debug('Received arguments: {}'.format(args))
    # Parsing of the config file to get primer3-settings
    config = configparser.ConfigParser()
    config.read(args.config.name)
    logging.info('Successfully read config {}'.format(args.config.name))
    if args.additionalFasta is not None:
        sequences = args.additionalFasta
        logging.info('Found additional file for primer generation {}'.format(
            sequences.name
        ))
    else:
        sequences = args.FastaFile
        logging.info(
            'No additional file with sequences for primer-generation specified')
    # extract sequence from sequences
    sequence, seq_id = parse_fasta(sequences, args.sequence)
    if not sequence:
        sys.exit("Could not find sequence with given ID-Prefix")
    # are there any settings for primer3?
    if config.has_section('primer3'):
        primer3_conf = config['primer3']
        logging.info('Found settings for primer3')
        logging.debug(
            'primer3-settings: {}'.format(' '.join(primer3_conf.keys())))
    else:
        primer3_conf = None
        logging.warning(
            'Config contains no settings for primer3, '
            'running primer3 with default settings.')

    # generate primers, dependent on sequence, primer3-configuration
    # and names of the files containing the found primers
    generate_primer(sequence, primer3_conf, args.primerfiles, config['default'])

    # is an already existing bowtie-index specified?
    if not args.index:
        # no index available, so we have to create our own one
        bowtie_index = setup_bowtie(args.FastaFile,
                                    args.loglevel == logging.DEBUG)
        logging.info("No existing index for bowtie specified")
    else:
        bowtie_index = args.index
        logging.info("Using existing bowtie-index")

    # any settings for bowtie?
    if config.has_section('bowtie'):
        bowtie_conf = config['bowtie']
        logging.info('Found settings for bowtie')
        logging.debug(
            'bowtie-settings: {}'.format(' '.join(bowtie_conf.keys())))
    else:
        bowtie_conf = None
        logging.info('Config contains no settings for bowtie')
    run_bowtie(bowtie_index, args.primerfiles, args.bowtie, bowtie_conf,
               args.loglevel == logging.WARNING)


def setup_bowtie(fasta_file, debug):
    """
    If no bowtie-index is specified we have to build it
    :type fasta_file: argparse.FileType
    """
    import re
    bowtie_index_dir = 'bowtie-index'
    # create new directory for the index,
    # no problem if specific folder already exists
    os.makedirs(bowtie_index_dir, exist_ok=True)
    # determine name for index from name of the
    # FASTA-file containing the sequences
    bowtie_index = "{index_dir}/{prefix}_bowtie".format(
        index_dir=bowtie_index_dir,
        prefix=
        re.split("/|\.", fasta_file.name)[-2])
    args = ["bowtie-build", fasta_file.name, bowtie_index]
    logging.info('bowtie-build command: {}'.format(args))
    if not debug:
        # We are not in debug-mode so no output will be shown
        subprocess.call(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        # debug-mode, show everything
        subprocess.call(args)
    return bowtie_index


def run_bowtie(bowtie_index, files_prefix, bowtie_exec, bowtie_config, silent):
    """
    Calls bowtie to execute the search for matches of the designed primers with
    other sequences.
    :param bowtie_config:
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

    # check whether settings for bowtie are available
    if bowtie_config is not None:
        args = parse_bowtie_config(args, bowtie_config)
    else:
        logging.warning('Config contains no settings for bowtie')

    logging.debug('Calling bowtie: {}'.format(args))

    if silent:
        res = subprocess.check_output(args, stderr=subprocess.PIPE).decode(
            'utf-8').split('\n')
    else:
        logging.info('Bowtie result summary:')
        res = subprocess.check_output(args).decode('utf-8').split('\n')

    logging.info('Bowtie result:\n{}'.format('\n'.join(res)))


def parse_bowtie_config(cmd, config):
    def is_int(x):
        """
        Checks whether passed value can be parsed as valid int
        :param x: Value to check
        :return: True if it can be parsed as int
        """
        try:
            return type(ast.literal_eval(x)) == int
        except SyntaxError:
            return False

    # TODO more settings; even custom
    for key in config.keys():
        # is the key defined as int?
        if is_int(config[key]):
            # are there any special restrictions defined?
            if key in bowtie_valid_setting.keys():
                (is_valid, mes) = bowtie_valid_setting[key]
                # test for special restrictions and show warning if they
                # do not hold
                if not is_valid(config.getint(key)):
                    logging.warning((
                        'Value for {conf} in config,'
                        ' resp. {bowtie} in bowtie: {val}'
                        ' is invalid. {msg}. Ignoring value.'
                    ).format(
                        conf=key, msg=mes,
                        bowtie=bowtie_options_translation[key],
                        val=config[key])
                    )
                    # skip this value since it is invalid
                    break
            # try to translate our settings name to a valid
            # option for bowtie
            if key in bowtie_options_translation.keys():
                cmd += [bowtie_options_translation[key],
                        config[key]]
            else:
                # option is not valid and will be ignored
                logging.warning(
                    ('Found unknown value in bowtie config: {}. '
                     'Ignoring value').format(key)
                )
        else:
            logging.warning(
                ('Found non-int value for {conf} in config,'
                 ' resp. {bowtie} in bowtie: {val}').format(
                    conf=key, bowtie=bowtie_options_translation[key],
                    val=config[key])
            )
    return cmd


def generate_primer(sequence, primer3_config, primer_file_prefix, config):
    """
    Calls the primer3-module with the settings and separates the results in
    left and right primer pairs.
    :param primer_file_prefix: prefix for the files where the primer pairs will
    be stored.
    :param sequence: Sequence template for the primers
    :param primer3_config: containing the settings for primer3
    :return: The separated primer-pairs.
    """
    primer3_config_dict = {}  # type: dict
    # are any specific settings for primer3?
    if primer3_config is not None:
        logging.info('Parsing primer3-settings from config')
        # if yes, we have to parse them
        for k in primer3_config.keys():
            try:
                value = ast.literal_eval(primer3_config[k])
                logging.debug('Evaluated {key} to {v}: {t}'.format(
                    key=k, v=value, t=type(value)
                ))
                primer3_config_dict.update({str(k).upper(): value})
            except SyntaxError:
                logging.warning(
                    ('Could not parse option {key} with value {v} '
                     'from primer3-config. '
                     'Ignoring value').format(key=k, v=primer3_config[k])
                )
    logging.debug('Final primer3-options-dictionary: {}'.format(
        [str(k) + ": " + str(primer3_config_dict[k]) for k in
         primer3_config_dict.keys()]
    ))  # set primer3-settings
    primer3.setP3Globals(primer3_config_dict)

    # remove any newlines or anything else like that
    sequence = sequence.replace('\n', '').replace('\r', '')
    # extract region inside sequence to generate primers
    seq_begin = config.getint('SEQUENCE_INCLUDED_BEGIN', -1)
    seq_end = config.getint('SEQUENCE_INCLUDED_END', -1)
    # is one of the values invalid?
    if seq_begin < 0 or seq_end < 0:
        sys.exit(
            ('ERROR:Found non-negative or invalid values for primer'
             ' generation concerning the region inside '
             'the sequence: (SEQUENCE_INCLUDED_BEGIN,SEQUENCE_INCLUDED_END)'
             ' = {}. If you want to generate primers for the '
             'whole sequence set both values to 0. Aborting').format(
                (config['SEQUENCE_INCLUDED_BEGIN'],
                 config['SEQUENCE_INCLUDED_END']))
        )

    # generate primers for the whole sequence?
    if seq_begin == seq_end == 0:
        logging.info(
            'Generating primers for the whole sequence as requested')
        res = primer3.bindings.designPrimers({
            'SEQUENCE_ID': 'mySequence',
            'SEQUENCE_TEMPLATE': sequence,
        })
    else:
        logging.info(
            'Generating primers for region {}'.format((seq_begin, seq_end)))
        res = primer3.bindings.designPrimers({
            'SEQUENCE_ID': 'mySequence',
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [seq_begin, seq_end]

        })

    primer_left = {}  # type: dict
    primer_right = {}  # type: dict
    for k in res.keys():
        line = k.split('_')
        # only applies to keys containing primer-sequences
        if len(line) >= 4:
            # only search for RIGHT or LEFT;
            # MIDDLE would be available as well but we only want the primer
            if line[1] in ['RIGHT'] and \
                    line[3] == 'SEQUENCE':
                primer_right.update({k: res[k]})
            elif line[1] in ['LEFT'] and \
                    line[3] == 'SEQUENCE':
                primer_left.update({k: res[k]})
    # write the found primer to their corresponding files
    logging.debug('Opening files to write primers')
    try:
        primerfile_left = open(
            '{prefix}_left.fas'.format(prefix=primer_file_prefix), 'x')
        primerfile_right = open(
            '{prefix}_right.fas'.format(prefix=primer_file_prefix), 'x')
    except FileExistsError:
        logging.error((
            'Found existing files with names ({prefix}_left.fas and/or '
            '{prefix}_right.fas) which should be used for the primerfiles. '
            'Please remove them and start the programm again').format(
            prefix=primer_file_prefix))
        sys.exit('Aborting')

    # create function for sorting the prefixes
    def primer_sort(x):
        return x.split('_')[2] + x.split('_')[1]

    logging.debug('Writing left primers to {}'.format(primerfile_left.name))
    for k in sorted(primer_left.keys(), key=primer_sort):
        # format in FASTA-style
        line = ">{}\n{}\n\n".format(k, primer_left[k])
        primerfile_left.write(line)
    primerfile_left.close()

    logging.debug(
        'Writing right primers to {}'.format(primerfile_right.name))
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
        ID of the sequence which should be given to primer3.
        If an additional FASTA-File for primer-generation is passed via '-a'
        it will be searched, otherwise the FastaFile.
        Prefix-matching is used for identification
        and the first hit will be used.
        """
    )
    parser.add_argument(
        "-a", "--additionalFasta", type=argparse.FileType('r'),
        default=None,
        help=
        """
        An additional file containing the sequence for which primer should be
        generated. First sequence inside of the file is taken if not specified
        otherwise via '-s'.
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
        '-d', '--debug',
        help="Print lots of debugging statements",
        action="store_const", dest="loglevel", const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-v', '--verbose',
        help="Be verbose",
        action="store_const", dest="loglevel", const=logging.INFO,
        default=logging.WARNING,
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
        # read one line to remove headerline from internal buffer
        seq_id_header = fasta_file.readline().split()[0][1:]
        logging.info(
            'No seq_id passed, taking first sequence from {} with id {}'.format(
                fasta_file.name, seq_id_header))
        for line in fasta_file:
            # read as long as no newline appears
            if line in ['\n', '\r\n']:
                break
            seq += line
    # a sequence-id has been passed to the programm
    else:
        logging.info('Partial sequence-id given: {}'.format(seq_id))
        for line in fasta_file:
            # begin of a new sequence and prefix-matching says true
            if line[0] == '>' and line.startswith(seq_id, 1):
                seq_id_header = line.split()[0][1:]
                logging.info(
                    'Found match of given sequence-id-prefix, '
                    'using sequence with id {}'.format(
                        seq_id_header))
                break
        for line in fasta_file:
            # again, read as long as no newline appears
            if line in ['\n', '\r\n']:
                break
            seq += line

    logging.debug('Successfully extracted sequence')
    return seq, seq_id_header


if __name__ == "__main__":
    main()
