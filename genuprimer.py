#!/usr/bin/env python3
"""
Author: tluettje@techfak.uni-bielefeld.de
License: AGPL v3
"""
import argparse
import ast
import configparser
import logging
import os
import re
import subprocess
import sys

import primer3

MANDATORY_VALUES = {'last_must_match', 'last_to_check', 'last_max_error'}
RESULT_HEADER = "FWD_ID,REV_ID,MATCH_ID,FWD,REV,START,STOP,LENGTH,EXP"
LOGGING_LEVEL = {'WARNING': logging.WARNING,
                 'INFO': logging.INFO,
                 'DEBUG': logging.DEBUG}

MAX_NUMBER_OF_MATCHES = 5

PRIMER3_OPTIONS = {}  # type: dict

BOWTIE_RUN_OPTIONS = {}  # type: dict

BOWTIE_PARSE_OPTIONS = {'LAST_MUST_MATCH': 3,
                        'LAST_TO_CHECK': 12,
                        'LAST_MAX_ERROR': 5,
                        'MAX_NUMBER_OF_MATCHES': 5}


def setup_logging(loglevel: str):
    """
    Setups the logging functionality, sets colors for different log levels
    and parses the loglevel chosen by user.
    :param loglevel: chosen loglevel by user
    :return: None
    """
    # setup logging and level as well as colors for logging
    # colors are only working on UNIX-Shells
    logging.addLevelName(logging.WARNING,
                         # format yellow
                         "\033[1;35m%s\033[1;0m" % logging.getLevelName(
                             logging.WARNING))
    logging.addLevelName(logging.ERROR,
                         # format red
                         "\033[1;31m%s\033[1;0m" % logging.getLevelName(
                             logging.ERROR))
    logging.addLevelName(logging.INFO,
                         # format green
                         "\033[1;32m%s\033[1;0m" % logging.getLevelName(
                             logging.INFO))
    logging.addLevelName(logging.DEBUG,
                         # format blue
                         "\033[1;36m%s\033[1;0m" % logging.getLevelName(
                             logging.DEBUG))
    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=LOGGING_LEVEL[loglevel])


def parse_config_and_parameters(args: argparse.Namespace,
                                config: configparser.ConfigParser):
    pass


def main():
    """
    Delegates all tasks to the other functions.
    """
    # parse all arguments, let the argparse-module do its wonderful work
    args = parse_arguments()
    # setup logging
    setup_logging(args.loglevel)
    logging.debug('Received arguments: {}'.format(args))
    # Parsing of the config file to get primer3-settings
    config = configparser.ConfigParser()
    config.read(args.config.name)
    logging.info('Successfully read config {}'.format(args.config.name))
    # check whether all mandatory values are present in config
    if not MANDATORY_VALUES.issubset(set(config['default'])):
        logging.error(
            'Not all mandatory values found in default-section '
            'of the config. Missing: {}'.format(
                list(map(lambda x: x.upper(),
                         set(MANDATORY_VALUES.difference(
                             set(config['default'])))))
            ))
        sys.exit(1)

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

    seq_included_region = extract_included_region(config['default'])

    """
    Existing primer pairs specified via -p/--primerfiles will be read and bowtie
    run against them.
    """
    if not args.keep_primer:
        logging.info('Generating new primer.')
        if args.additionalFasta is not None:
            sequences = args.additionalFasta
            logging.info(
                'Found additional file for primer generation {}'.format(
                    sequences.name
                ))
            logging.warning(
                'It is not possible to say whether a match found '
                'by bowtie is expected or not if an additional file for '
                'primer generation is passed.')
        else:
            sequences = args.FastaFile
            logging.info(
                'No additional file with sequences for primer-generation '
                'specified')
        """
        Extract sequence from FASTA file and get exact id of it as well
        """
        sequence, seq_id = parse_fasta(sequences, args.sequence)
        if not sequence:
            logging.error(
                "Could not find sequence with given ID-Prefix. Aborting")
            sys.exit(1)
        logging.info('Successfully extracted sequence')

        """
        Generate primer pairs, depending on extracted sequence, primer3
        configuration and specified region for which primer shall be generated.
        """
        primer_dict = generate_primer(sequence, primer3_conf, args.primerfiles,
                                      seq_included_region)
    else:
        logging.info('Trying to parse existing primer from files specified'
                     ' via -p/--primerfiles')
        logging.warning('Calculations whether a hit is expected or not '
                        'are based on the settings which would be used for '
                        'normal primer generation: SEQUENCE_INCLUDED_* for '
                        'the region, PRIMER_PRODUCT_SIZE_RANGE for size of '
                        'a possible insert and -s/--sequence for prefix-'
                        'matching with the reported id by bowtie in which '
                        'sequence the hit was found.')
        # store given sequence-parameter for prefix-matching with bowtie-results
        if args.sequence is not None:
            seq_id = args.sequence
        else:
            seq_id = None

        primer_dict = parse_existing_primer(args.primerfiles)
        # we did not generate any primer, therefore we do not have any sequence
        # id to check whether a match is expected or not
    # is an already existing bowtie-index specified?
    if not args.index:
        # no index available, so we have to create our own one
        bowtie_index = setup_bowtie(args.FastaFile,
                                    args.loglevel == logging.DEBUG, args.bowtie)
        logging.info("No existing index for bowtie specified")
    else:
        bowtie_index = args.index
        logging.info("Using existing bowtie-index")

    bowtie_result = run_bowtie(bowtie_index, args.primerfiles, args.bowtie,
                               args.loglevel == logging.WARNING,
                               product_size_range, args.bowtie_output)

    # create empty list for results
    results = {}
    for primer_tuple in bowtie_result:
        current_key, res = parse_bowtie_result(primer_tuple, config['default'],
                                               primer_dict, seq_included_region,
                                               args.additionalFasta is not None,
                                               seq_id, args.keep_primer)
        if res is not None:
            results.setdefault(current_key, []).append(res)
        else:
            continue

    # get value defined in config or fallback otherwise
    max_number_matches = config['default'].getint('MAX_NUMBER_OF_MATCHES',
                                                  MAX_NUMBER_OF_MATCHES)
    logging.info('Maximal number of allowed matches set to {}'.format(
        max_number_matches))

    # set output either to STDOUT or output file if specified by user
    if args.output is None:
        logging.info('No file for output specified, writing to STDOUT.')
        output = sys.stdout
    else:
        logging.info('Output file specified, writing results to {}'.format(
            args.output.name))
        output = args.output

    # store intermediate all results which would be printed in output
    printable_res = []
    for key in sorted(results.keys()):
        matches = results[key]
        if len(matches) > max_number_matches:
            logging.debug(
                'Not printing results for {} because it has {} matches'.format(
                    key, len(matches)
                ))
        else:
            # add to intermediate results
            printable_res.append(matches)

    output.write(RESULT_HEADER + '\n')
    for matches in sorted(printable_res, key=len):
        # write results
        output.write('\n'.join(matches))
        # final newline at end of results
        output.write('\n')


def parse_existing_primer(prefix: str) -> dict:
    """
    This function reads the files containing the custom primer pairs and stores
    them in the same way newly generated would be stored. Only complete pairs
    are stored. If one file contains m and the other n sequences, given n < m
    only n pairs would be returned and the m-n discarded.
    :param prefix: Prefix for files containing the primer, expected to be
    readable from prefix_{left,right}.fas
    :return: Dictionary containing all primer, accessible via sorted
    concatenation of their names
    """
    left_name = "{}_left.fas".format(prefix)
    right_name = "{}_right.fas".format(prefix)
    try:
        left_primer = open(left_name, 'r')
        right_primer = open(right_name, 'r')
    except FileNotFoundError:
        logging.error(
            'Could not find or read one or both of the files '
            'containing the primers that should be used for this run '
            'instead of generating new one.\n'
            'Looking for files: {} and {}. Aborting'.format(
                left_name, right_name
            ))
        sys.exit(1)

    primer_left_lines = left_primer.readlines()
    primer_right_lines = right_primer.readlines()
    left_primer.close()
    right_primer.close()

    primer_left_list = []
    sequence = seq_id = ""
    for left_line in primer_left_lines:
        # begin of new sequence
        if left_line[0] == '>':
            if sequence != "":
                # we were previously reading a sequence therefore store it
                # before starting with a new one

                # sanitize sequence from newlines
                sequence = sequence.replace('\n', '').replace('\r', '')
                primer_left_list.append((seq_id, sequence))
                sequence = ""
            # extract id from header
            seq_id = left_line.split('>')[1].split()[0]
        else:
            sequence += left_line
    # append final sequence to list
    if sequence != "":
        sequence = sequence.replace('\n', '').replace('\r', '')
        primer_left_list.append((seq_id, sequence))

    primer_right_list = []
    sequence = seq_id = ""
    for right_line in primer_right_lines:
        # begin of new sequence
        if right_line[0] == '>':
            if sequence != "":
                # we were previously reading a sequence therefore store it
                # before starting with a new one
                sequence = sequence.replace('\n', '').replace('\r', '')
                # sanitize sequence from newlines
                primer_right_list.append((seq_id, sequence))
                sequence = ""
            # extract id from header
            seq_id = right_line.split('>')[1].split()[0]
        else:
            sequence += right_line
    if sequence != "":
        sequence = sequence.replace('\n', '').replace('\r', '')
        primer_right_list.append((seq_id, sequence))

    logging.debug('Extracted following primer from {}: {}'.format(
        left_name, ' ,'.join(map(str, primer_left_list))
    ))
    logging.debug('Extracted following primer from {}: {}'.format(
        right_name, ' ,'.join(map(str, primer_right_list))
    ))
    primer_dict = {}  # type: dict
    for (l_id, l_seq), (r_id, r_seq) in zip(primer_left_list,
                                            primer_right_list):
        l_seq = l_seq.replace('\n', '').replace('\r', '')
        r_seq = r_seq.replace('\n', '').replace('\r', '')
        primer_dict.update({tuple(sorted((l_id, r_id))): (l_seq, r_seq)})
    return primer_dict


def extract_fasta_seq(fasta_file: '_io.TextIOWrapper') -> str:
    """
    Gets a file buffer pointing at the beginning of a sequence. The file is read
    and the line content added to an internal buffer until a the sequence ends
    which is either indicated by newline or the beginning of a new sequence
    :param fasta_file: IO-wrapper of the file
    :return: extracted sequence sanitized from newlines
    """
    seq = ""
    for line in fasta_file:
        # read as long as no newline appears
        if line in ['\n', '\r\n'] or line[0] == '>':
            break
        seq += line
    return seq.replace('\n', '').replace('\r', '')


def extract_included_region(config: configparser.SectionProxy) -> tuple:
    """
    This functions extracts the region for which primer shall be generated
    from the config file.
    :param config: contains the values.
    :return: extracted values formatted as int
    """
    if config is not None:
        seq_begin = config.getint('SEQUENCE_INCLUDED_BEGIN', -1)
        seq_end = config.getint('SEQUENCE_INCLUDED_END', -1)
    else:
        seq_begin = 0
        seq_end = 0
    # is one of the values invalid?
    if seq_begin < 0 or seq_end < 0:
        logging.error(
            ('Found negative or invalid values for the region'
             'for which the primer will be generated or have been generated if '
             'custom primers are used via --keep-primer.'
             '(SEQUENCE_INCLUDED_BEGIN,SEQUENCE_INCLUDED_END)'
             ' = {}. If you want to generate primers for the '
             'whole sequence set both values to 0. Aborting').format(
                (config['SEQUENCE_INCLUDED_BEGIN'],
                 config['SEQUENCE_INCLUDED_END']))
        )
        sys.exit(1)
    return seq_begin, seq_end


def parse_bowtie_result(primer_tuple: tuple,
                        error_values: configparser.SectionProxy,
                        primer_dict: dict,
                        seq_included_region: tuple, additional_fasta: bool,
                        seq_id: str,
                        keep_primer: bool) -> tuple:
    """
    Parses the bowtie result and presents them in a csv-style containing
    the most important information.
    :param primer_tuple: Current tuple (forward and reverse) which bowtie
    results are now checked.
    :param error_values: Dictionary containing the values which are necessary
    to check for valid matches.
    :param primer_dict: Dictionary containing all primer pairs, accessible via
    sorted concatenation of fwd-primer-id and rev-primer-id.
    :param seq_included_region: Start and stop of the region for which the
    primer have been generated.
    :param additional_fasta: Whether the extracted sequence's origin is an
    additional file not indexed by bowtie.
    :param seq_id: The whole id of the sequence for which primer have been
    generated or a prefix defined by user if custom primer are used.
    :param keep_primer: Whether custom primer pairs were read from files and
    no new ones were generated.
    :param seq_included_region: Tuple of the specified region for which primer
    have been produced; used to check whether a match reported by bowtie is
    expected or not
    :return tuple of key and result for this key, where the key is the same from
    primer_dict
    """
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

    def is_significant(values):
        # if last entry is a number it represents the number of matches at
        # the end of the primer; must not be lower than given value
        if values[-1].isdigit():
            if int(values[-1]) < last_must_match:
                logging.debug(
                    'Failed because last_must_match not fulfilled: {}'.format(
                        values))
                return False
            # only one number at the end and it is greate than last_must_match
            if len(values) == 1:
                return True
        # if last entry is a char there is no match, unless we do not care about
        # the entries at the end
        elif values[-1].isalpha() and last_must_match != 0:
            logging.debug(
                'Failed because last char is str, '
                'therefore last n bases cannot match: {}'.format(values))
            return False

        logging.debug('Checking non-trivial: {}'.format(values))
        # calculate the last entries to check
        # index to iterate over result values
        i = 1
        # counted number of substitutions/mismatches
        number_of_subs = 0
        # how many bases, when expanding numeral values are already processed
        bases_processed = 0
        # can we stop checking?
        while bases_processed <= last_to_check:
            current = values[-i]
            if current.isdigit():
                bases_processed += int(current)
            else:
                number_of_subs += 1
                bases_processed += 1
            i += 1
        if number_of_subs < last_max_error:
            return True
        else:
            logging.debug(
                'Failed because of {num} subs in last {last} '
                'bases, {max} allowed'.format(
                    num=number_of_subs, last=last_to_check, max=last_max_error
                ))
        return number_of_subs < last_max_error

    # if bowtie result line is so short it means no hit has been found
    # and obviously a non existent match cannot be significant
    if len(primer_tuple[0].split()) <= 12 or len(primer_tuple[1].split()) <= 12:
        # if bowtie result line is so short it means no hit has been found
        # and obviously a non existent match cannot be significant
        return None, None
    # split bowtie results for forward and reverse primer of one pair
    [left_split, right_split] = list(map(str.split, primer_tuple))
    # get string representation of mismatch bases
    left_match, right_match = left_split[12], right_split[12]
    # get id of each primer pair
    left_name, right_name = left_split[0], right_split[0]
    # split numerical and alphabetic values
    left_res = re.split('(\d+)', left_match.split(':')[2])
    right_res = re.split('(\d+)', right_match.split(':')[2])

    def remove_empty_mismatch(x):
        """
        Used to remove any non important information from the string match
        representation of bowtie. Empty ones are not informative and the usage
        of 0 does not make any sense me since it means: 0 consecutive bases
        matched.
        :param x: A digit or char from the match string representation
        :return: Boolean value whether it should be removed or not
        """
        return x not in ['', '0', 0]

    # sanitize string representations from empty values
    left_res = list(filter(remove_empty_mismatch, left_res))
    right_res = list(filter(remove_empty_mismatch, right_res))
    # check whether the hit reported bowtie is significant according to users
    # preferences
    if is_significant(left_res) and is_significant(right_res):
        # extract information from one of the results
        infos = left_split
        # generate key for dictionary containing all primer
        # for current bowtie results
        current_dict_key = tuple(sorted((left_name, right_name)))
        # get current sequences
        current_pair = primer_dict.get(current_dict_key)

        if additional_fasta:
            """
            If primer have been generated using an additional file we cannot
            compare positional information with the fasta file used in bowtie
            therefore every hit will be not expected.
            """
            expected_hit = False
        elif keep_primer:
            """
            Primer have not been generated by us but user has been warned that
            calculations whether a hit is expected or will be based on given
            values nevertheless. INCLUDED_REGION_* values are used.
            But in contrast to normal computation we cannot expect that the id
            of a sequence where bowtie has found a match will be same as the one
            from our primers since we only have the value from the user.
            We therefore look whether the given id is a prefix of the reported
            id.
            """
            expected_hit = seq_included_region[0] <= int(infos[3]) <= int(
                infos[7]) <= seq_included_region[1] and infos[2].startswith(
                seq_id)
        else:
            """
            A match is expected if its start and end position are inside
            SEQUENCE_INCLUDED_REGION used for primer generation and if
            the id of the sequence where the match was found is the
            same as the one for which the primer were generated.
            """
            expected_hit = seq_included_region[0] <= int(infos[3]) <= int(
                infos[7]) <= seq_included_region[1] and infos[2] == seq_id

        # format results to result format
        res = ('{fwd}{sep}{rev}{sep}{gi}{sep}{left_primer}{sep}{right_primer}'
               '{sep}{start}{sep}{stop}{sep}{size}{sep}{exp}').format(
            fwd=left_name, rev=right_name, sep=',', gi=infos[2],
            left_primer=current_pair[0],
            right_primer=current_pair[1], start=infos[3],
            # bowtie sets stop to the first of the rev primer, but we want to
            # the length of the whole region enclosed by the primer including
            # themselves
            stop=int(infos[7]) + len(current_pair[1]),
            size=infos[8], exp=1 if expected_hit else 0
        )
        return current_dict_key, res
    else:
        return None, None


def setup_bowtie(fasta_file: '_io.TextIOWrapper', debug: bool,
                 bowtie_exec: str) -> str:
    """
    If no bowtie-index is specified we have to build it.
    :param fasta_file: io-wrapper of the file with sequences that shall be
    indexed, needed to determine name of index
    :param debug: whether debug logging is set on or off.
    :param bowtie_exec: str containing the path to bowtie executable.
    bowtie-build is supposed to be in the same folder.
    :return location of the build index
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
        prefix=re.split("/|\.", fasta_file.name)[-2])
    bowtie_build = bowtie_exec + '-build'
    args = [bowtie_build, fasta_file.name, bowtie_index]
    logging.info('bowtie-build command: {}'.format(args))
    if not debug:
        # We are not in debug-mode so no output will be shown
        subprocess.call(args, stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL)
    else:
        # debug-mode, show everything
        subprocess.call(args)
    return bowtie_index


def run_bowtie(bowtie_index: str, files_prefix: str, bowtie_exec: str,
               silent: bool, size_range: tuple,
               bowtie_output: bool) -> list:
    """
    Calls bowtie to execute the search for matches of the designed primers with
    other sequences.
    :param bowtie_index: location of the index for bowtie
    :param files_prefix: the prefix of the files containing the primer
    :param bowtie_exec: bowtie-executable, either default or defined by user
    :param silent: whether we run in silent, decides whether bowtie output
    shall be shown
    :param size_range: given size range of the primer, we therefore only look
    for inserts of this size range
    :param bowtie_output: whether output of bowtie shall be written to STDERR
    :return: list of tuples consisting of the hits
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
            if type(tmp) == list and isinstance(tmp[0], int) and isinstance(
                tmp[1], int):
                args += ['--minins', str(tmp[0]), '--maxins', str(tmp[1])]

        except (ValueError, SyntaxError):
            # Same error has already been reported during primer generation
            pass

        logging.debug(
            'Product size range was specified for primer3'
            '. Applying values for bowtie insert size: {}'.format(size_range)
        )
    logging.info('Calling bowtie: {}'.format(args))
    if silent:
        args += ['--quiet']
        res = subprocess.check_output(args).decode('utf-8').split('\n')
    else:
        logging.info('Bowtie result summary:')
        res = subprocess.check_output(args).decode('utf-8').split('\n')

    if bowtie_output:
        logging.info('Printing bowtie result to STDERR as requested by '
                     '--show-bowtie')
        print('\n'.join(res), file=sys.stderr)
        sys.stderr.flush()

    # Each match is described in two lines since FWD and REV have to match.
    res_tuple = [
        (res[i], res[i + 1]) for i in range(0, len(res) - 1, 2)]
    return res_tuple


def generate_primer(sequence: str, primer3_config: configparser.SectionProxy,
                    primer_file_prefix: str,
                    seq_included_region: tuple) -> dict:
    """
    Calls the primer3-module with the settings and separates the results in
    left and right primer pairs.
    :param sequence: Sequence template for the primers
    :param primer3_config: containing the settings for primer3
    :param primer_file_prefix: prefix for the files where the primer pairs will
    be stored.
    :param seq_included_region: tuple containing positions for which
    primer shall be generated
    :return: A dictionary containing all primer pairs and their names
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

    (seq_begin, seq_end) = seq_included_region
    # generate primers for the whole sequence?
    if seq_begin == seq_end == 0:
        logging.warning(
            'Generating primers for the whole sequence as '
            'requested since SEQUENCE_INCLUDED_END and _BEGIN are set to 0')
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
            # give start of sequence and length
            'SEQUENCE_INCLUDED_REGION': [seq_begin, seq_end - seq_begin]

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

    def extract_number(x: str) -> int:
        """
        Extracts the number of the name of a generated primer
        :param x: primer-id
        :return: included number
        """
        return int(x.split('_')[2])

    primer_dict = {}

    logging.debug(
        'Writing primers to files: {}, {}'.format(primerfile_left.name,
                                                  primerfile_right.name))
    for left_key, right_key in zip(
        sorted(primer_left.keys(), key=extract_number),
        sorted(primer_right.keys(), key=extract_number)
    ):
        # format in FASTA-style
        left_line = ">{}\n{}\n\n".format(left_key, primer_left[left_key])
        right_line = ">{}\n{}\n\n".format(right_key, primer_right[right_key])
        primerfile_left.write(left_line)
        primerfile_right.write(right_line)
        primer_dict.update({tuple(sorted((left_key, right_key))): (
            primer_left[left_key], primer_right[right_key])})

    primerfile_left.close()
    primerfile_right.close()

    return primer_dict


def parse_arguments() -> argparse.Namespace:
    """
    This function parses the commandline arguments via the argparse-module,
    which additionally generates a help-hook, if a parameter is passed wrong.
    :return: A parser-object containing all parsed values.
    """
    parser = argparse.ArgumentParser(
        description=
        """
        genuprimer is able to generate new primer by calling primer3 and
        afterwards validate their uniqueness among other sequences by calling
        bowtie. The same can be accomplished for already existing primer pairs.
        """
        , formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
        Partial ID of the sequence for which the primer shall be or have been
        generated.
        Used for primer generation to extract the chosen sequence from FastaFile
        and during validation by bowtie if --keep-primer has been set.
        If this option is skipped the first sequence from FastaFile or
        additionalFasta, if omitted, is used.
        """
    )
    parser.add_argument(
        "-a", "--additionalFasta", type=argparse.FileType('r'),
        default=None,
        help=
        """
        An additional file containing the sequence for which primer shall be
        generated. First sequence inside of the file is taken if not specified
        otherwise via '-s'.
        """
    )
    parser.add_argument(
        "-p", "--primerfiles", type=str, default="primer3", help=
        """
        Prefix for the files where the primer pairs will be stored if new ones
        are generated. The suffixes will be '_left.fas' and '_right.fas'.
        Therefore if the default remains, their names will be
        'primer3_left.fas' and 'primer3_right.fas'. If using --keep-primer
        the custom files must follow this convention.
        """
    )
    parser.add_argument(
        '-i', '--index', type=str, default=None, help=
        """
        Specify an existing bowtie-index otherwise a new one will be generated
        for FastaFile. This option is directly forwarded to bowtie.
        """
    )
    parser.add_argument(
        '-v', '--verbose',
        help="Be verbose by showing INFO messages.",
        action="store_const", dest="loglevel", const='INFO',
        default='WARNING',
    )
    parser.add_argument(
        '-d', '--debug',
        help="Print lots of DEBUG messages.",
        action="store_const", dest="loglevel", const='DEBUG',
        default='WARNING',
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
        default=None, help=
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
        The bowtie executable if not in PATH. If needed bowtie-build is expected
        to be found via appending '-build' to bowtie.
        """
    )
    parser.add_argument(
        "--show-bowtie", dest='bowtie_output', action='store_true',
        help=
        """
        Set this option to show the original results of bowtie, written to
        STDERR.
        """
    )
    parser.add_argument(
        "--keep-primer", dest='keep_primer', action='store_true',
        help=
        """
        Set this option to start another run with the same primers from last run
        or some custom one.
        """
    )
    parser.set_defaults(keep_primer=False, bowtie_output=False)
    return parser.parse_args()


def parse_fasta(fasta_file: '_io.TextIOWrapper', seq_id: str) -> tuple:
    """
    Parses the submitted FASTA-File and extracts the sequence for primer3.
    :param seq_id: identifier of the sequence, for which primers shall
    be generated
    :param fasta_file: already readable-opened file which
    contains all the sequences
    """
    seq_id_header = ""
    # no sequence-id specified, therefore the first one is taken
    if not seq_id:
        # read one line to remove headerline from internal buffer
        seq_id_header = fasta_file.readline().split()[0][1:]
        logging.info(
            'No seq_id passed, taking first sequence from {} with id {}'.format(
                fasta_file.name, seq_id_header))
        seq = extract_fasta_seq(fasta_file)
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
        seq = extract_fasta_seq(fasta_file)

    return seq, seq_id_header


if __name__ == "__main__":
    main()
