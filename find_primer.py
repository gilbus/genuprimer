import argparse
import ast
import configparser
import logging
import os
import subprocess
import sys
from argparse import RawTextHelpFormatter

import primer3

MANDATORY_VALUES = {'last_must_match', 'last_to_check', 'last_max_error'}
RESULT_HEADER = "PAIR_NR,ID,FWD,REV,START,STOP,LENGTH"


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
    logging.addLevelName(logging.ERROR,
                         # format red
                         "\033[1;31m%s\033[1;0m" % logging.getLevelName(
                             logging.ERROR))
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
    if args.config is not None:
        config.read(args.config.name)
        logging.info('Successfully read config {}'.format(args.config.name))
    else:
        logging.warning('No config file given. Using default values.')
    if not MANDATORY_VALUES.issubset(set(config['default'])):
        logging.error('Not all mandatory values found in default-section '
                      'of the config. Missing: {}'.format(
            list(map(lambda x: x.upper(),
                     set(MANDATORY_VALUES.difference(set(config['default'])))))
        ))
        sys.exit(1)
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
        logging.error("Could not find sequence with given ID-Prefix. Aborting")
        sys.exit(1)
    logging.debug('Successfully extracted sequence')

    # length of the insert between the primer pair
    product_size_range = None
    # are there any settings for primer3?
    if config.has_section('primer3'):
        primer3_conf = config['primer3']
        logging.info('Found settings for primer3')
        logging.debug(
            'primer3-settings: {}'.format(' '.join(primer3_conf.keys())))
        if 'primer_product_size_range' in primer3_conf.keys():
            product_size_range = primer3_conf['primer_product_size_range']
    else:
        primer3_conf = None
        logging.warning(
            'Config contains no settings for primer3, '
            'running primer3 with default settings.')

    if args.config is None:
        general_conf = None
    else:
        general_conf = config['default']
    # generate primers, dependent on sequence, primer3-configuration
    # and names of the files containing the found primers
    primer_pairs_list = generate_primer(sequence, primer3_conf,
                                        args.primerfiles, general_conf)

    # is an already existing bowtie-index specified?
    if not args.index:
        # no index available, so we have to create our own one
        bowtie_index = setup_bowtie(args.FastaFile,
                                    args.loglevel == logging.DEBUG)
        logging.info("No existing index for bowtie specified")
    else:
        bowtie_index = args.index
        logging.info("Using existing bowtie-index")

    bowtie_result = run_bowtie(bowtie_index, args.primerfiles, args.bowtie,
                               args.loglevel == logging.WARNING,
                               product_size_range)

    results = [RESULT_HEADER]
    for primer_tuple in bowtie_result:
        res = parse_bowtie_result(primer_tuple, config['default'],
                                  primer_pairs_list)
        results.append(res)

    if args.output.name == 'STDOUT':
        logging.info('No file for output specified, writing to STDOUT.')
        print('\n'.join(results))
    else:
        logging.info('Output file specified, writing results to {}'.format(
            args.output.name))
        for line in results:
            args.output.write('{}\n'.format(line))


def parse_bowtie_result(tuple, error_values, primer_pairs_list):
    try:
        last_must_match = int(error_values['last_must_match'])
        last_to_check = int(error_values['last_to_check'])
        last_max_error = int(error_values['last_max_error'])
    except ValueError:
        logging.error(
            ('One settings of {} from the config '
             'could not be parsed as integer. Aborting').format(
                ', '.join(MANDATORY_VALUES)
            ))
        sys.exit(1)

    def is_valid(values):
        logging.debug('checking: {}'.format(values))
        # if last entry is a number it represents the number of matches at
        # the end of the primer; must not be lower than given value
        if type(values[-1]) == int and int(values[-1]) < last_must_match:
            logging.debug('Failed because last_must_match not fulfilled')
            return False
        # if last entry is a char there is no match, unless we do not care about
        # the entries at the end
        elif values[-1].isalpha() and last_must_match != 0:
            logging.debug('Failed because last char is str')
            return False

        logging.debug(
            'last {} bases match for {}'.format(last_must_match, values))
        # calculate the last entries to check
        i = 1
        number_of_subs = 0
        while i <= last_to_check:
            logging.debug(i)
            current = values[-i]
            if current.isdigit():
                i += int(current)
            else:
                number_of_subs += 1
            i += 1
        logging.debug(
            'Counted {} errors in {}'.format(number_of_subs, values))
        return number_of_subs < last_max_error

    left_match, right_match = tuple[0].split()[12], tuple[1].split()[12]
    left, right = left_match.split(':')[2:], right_match.split(':')[2:]
    if is_valid(left) and is_valid(right):
        infos = tuple[0].split('\t')
        pair_number = int(infos[0].split('_')[2])
        current_pair = primer_pairs_list[pair_number]
        res = ('PRIMER_{num}{sep}{gi}{sep}{left_primer}{sep}{right_primer}'
               '{sep}{start}{sep}{stop}{sep}{size}').format(
            num=pair_number, sep=',', gi=infos[2], left_primer=current_pair[0],
            right_primer=current_pair[1], start=infos[3], stop=infos[7],
            size=infos[8]
        )
        return res


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


def run_bowtie(bowtie_index, files_prefix, bowtie_exec, silent, size_range):
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
    if size_range is not None:
        # value is specified, try to parse it as list
        try:
            tmp = ast.literal_eval(size_range)
            if type(tmp) == list and type(tmp[0]) == type(tmp[1]) == int:
                args += ['--minins', str(tmp[0]), '--maxins', str(tmp[1])]

        except (ValueError, SyntaxError):
            # Same error has already been reported during primer generation
            pass

        logging.debug('Product size range was specified for primer3'
                      '. Applying values for bowtie insert size: {}'.format(
            size_range
        ))
    logging.debug('Calling bowtie: {}'.format(args))
    if silent:
        res = subprocess.check_output(args, stderr=subprocess.PIPE).decode(
            'utf-8').split('\n')
    else:
        logging.info('Bowtie result summary:')
        res = subprocess.check_output(args).decode('utf-8').split('\n')

    logging.info('Bowtie result:\n{}'.format('\n'.join(res)))
    res_tuple = [
        (res[i], res[i + 1]) for i in range(0, len(res) - 1, 2)]
    logging.debug('Result tuple from bowtie: {}'.format(res_tuple))
    return res_tuple


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
            except (SyntaxError, ValueError):
                logging.warning(
                    ('Could not parse option {key} with value {v} '
                     'from primer3-config. '
                     'Ignoring value').format(key=k, v=primer3_config[k])
                )
    logging.debug('Final primer3-options-dictionary: {}'.format(
        [str(k) + ": " + str(primer3_config_dict[k]) for k in
         primer3_config_dict.keys()]
    ))  # set primer3-settings
    try:
        primer3.setP3Globals(primer3_config_dict)
    except TypeError as message:
        logging.error(
            'Settings for primer3 contain following error: {}\nAborting'.format(
                message
            ))
        sys.exit(1)

    # remove any newlines or anything else like that
    sequence = sequence.replace('\n', '').replace('\r', '')
    # extract region inside sequence to generate primers if given inside config
    if config is not None:
        seq_begin = config.getint('SEQUENCE_INCLUDED_BEGIN', -1)
        seq_end = config.getint('SEQUENCE_INCLUDED_END', -1)
    else:
        seq_begin = 0
        seq_end = 0
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
            'Generating primers for the whole sequence as requested or by default')
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
    primerfile_left = open(
        '{prefix}_left.fas'.format(prefix=primer_file_prefix), 'w')
    primerfile_right = open(
        '{prefix}_right.fas'.format(prefix=primer_file_prefix), 'w')

    # create function for sorting the prefixes
    def primer_sort(x):
        return x.split('_')[2] + x.split('_')[1]

    primer_pairs_list = []
    logging.debug('Writing left primers to {}'.format(primerfile_left.name))
    for k in sorted(primer_left.keys(), key=primer_sort):
        primer_pairs_list.append((primer_left[k],))
        # format in FASTA-style
        line = ">{}\n{}\n\n".format(k, primer_left[k])
        primerfile_left.write(line)
    primerfile_left.close()

    logging.debug(
        'Writing right primers to {}'.format(primerfile_right.name))
    for n, k in enumerate(sorted(primer_right.keys(), key=primer_sort)):
        primer_pairs_list[n] = (primer_pairs_list[n][0], primer_right[k])
        # format in FASTA-style
        line = ">{}\n{}\n\n".format(k, primer_right[k])
        primerfile_right.write(line)
    primerfile_right.close()

    return primer_pairs_list


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
        default=None, help=
        """
        Configfile with various parameters which are
        passed through to primer3. Has to include a
        'default'-section
        at top of the file, a 'primer3' section with the settings and a
        'bowtie' section.
        A default config is distributed with this programm.
        Non-working example:
        [default]
        SEQUENCE_INCLUDED_BEGIN = 30
        # other settings

        [bowtie]
        MaxInsSize = 100
        # more settings

        [primer3]
        PRIMER_OPT_SIZE = 14
        # more settings
        """, required=True)
    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'),
        default='STDOUT', help=
        """
        File where the results should be stored. Default is printing to STDOUT.
        Results are written as comma separated values (.csv).
        """
    )
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
    seq_id_header = ""
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

    return seq, seq_id_header


if __name__ == "__main__":
    main()
