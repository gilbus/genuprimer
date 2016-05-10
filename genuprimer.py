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

# Constants
RESULT_HEADER = "FWD_ID,REV_ID,MATCH_ID,FWD,REV,START,STOP,LENGTH,EXP"
LOGGING_LEVEL = {'WARNING': logging.WARNING,
                 'INFO': logging.INFO,
                 'DEBUG': logging.DEBUG}
CONFIG_REGION_KEYS = {'TARGET_POSITION_BEGIN': None,
                      'TARGET_POSITION_END': None,
                      'PRIMER_PRODUCT_SIZE_MIN': None,
                      'PRIMER_PRODUCT_SIZE_MAX': None}

primer3_insert_size = ()  # type: tuple

primer3_insert_pos = ()  # type: tuple

bowtie_parse_options = {'LAST_MUST_MATCH': 3,
                        'LAST_TO_CHECK': 12,
                        'LAST_MAX_ERROR': 5,
                        'LIMIT_NUMBER_OF_MATCHES': 5}

primer3_options = {}  # type: dict

# various runtime parameters and their default values
runtime_parameters = {'fasta_file': None,
                      'seq_id': '',
                      'config': 'genuprimer.conf',
                      'additional_fasta': None,
                      'index': 'bowtie-index',
                      'output': sys.stdout,
                      'prefix': 'genuprimer',
                      'bowtie': 'bowtie',
                      'keep_primer': False,
                      'show_bowtie_output': False}


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
    # convert args namespace to normal dictionary
    args = vars(args)
    if config:
        logging.info('Parsing configfile')
        if config.has_section('default'):
            logging.info('Found [default]-config in configfile')
            section = config['default']  # type: configparser.SectionProxy
            for key in section:  # type: str
                if key.upper() in bowtie_parse_options.keys():
                    logging.debug(('Found key {k} with value {v} in '
                                   '[default]-config').format(
                        k=key.upper(), v=section[key])
                    )
                    # get key from config and if it fails during conversion
                    # to an integer use the default value
                    try:
                        value = section.getint(
                            key.upper(),
                            fallback=bowtie_parse_options[key.upper()])
                    except ValueError:
                        # value is not a number, use fallback value
                        value = bowtie_parse_options[key.upper()]
                        logging.warning(
                            'Key {k} in [default]-config is not a number'.format(
                                k=key.upper()
                            ))
                    bowtie_parse_options[key.upper()] = value
                    logging.debug('New value for {k}: {v}'.format(
                        k=key.upper(), v=value))

                if key.upper() in CONFIG_REGION_KEYS.keys():
                    logging.debug(('Found key {k} with value {v} in '
                                   '[default]-config').format(
                        k=key.upper(), v=section[key])
                    )
                    # get key from config and if it fails during conversion
                    # to an integer use the default value
                    try:
                        value = section.getint(key.upper(), -1)
                    except ValueError:
                        value = CONFIG_REGION_KEYS[key.upper()]
                        logging.warning(
                            'Key {k} in [default]-config is not a number'.format(
                                k=key.upper()
                            ))
                    CONFIG_REGION_KEYS[key.upper()] = value
                    logging.debug('New value for {k}: {v}'.format(
                        k=key.upper(), v=value
                    ))
        if config.has_section('primer3'):
            logging.info('Found [primer3]-config in configfile')
            section = config['primer3']
            # if yes, we have to parse them
            for k in section.keys():
                try:
                    value = ast.literal_eval(section[k])
                    logging.debug('Evaluated {key} to {v}: {t}'.format(
                        key=k, v=value, t=type(value)
                    ))
                    primer3_options.update({str(k).upper(): value})
                except (SyntaxError, ValueError):
                    logging.warning(
                        ('Could not parse option {key} with value {v} '
                         'from [primer3] section.'
                         'Ignoring value').format(key=k, v=section[k])
                    )
        logging.info('Finished parsing of configfile')
        logging.debug('Current config values: {}'.format(
            str(CONFIG_REGION_KEYS) + str(bowtie_parse_options)
        ))

    logging.info('Parsing commandline arguments')
    runtime_parameters['show_bowtie_output'] = args['show_bowtie_output']
    runtime_parameters['keep_primer'] = args['keep_primer']
    for key in bowtie_parse_options:
        if args[key]:
            bowtie_parse_options[key] = args[key]
            logging.debug('New value for {k}: {v}'.format(
                k=key, v=bowtie_parse_options[key]
            ))
    if args['size']:
        global primer3_insert_size
        primer3_insert_size = tuple(args['size'])
        logging.debug('New value for PRIMER_PRODUCT_SIZE: {v}'.format(
            v=primer3_insert_size
        ))
    if args['pos']:
        global primer3_insert_pos
        primer3_insert_pos = tuple(args['pos'])
        logging.debug('New value for PRIMER_INSERT_POSITION: {v}'.format(
            v=primer3_insert_pos
        ))
    for key in runtime_parameters:
        if key in args.keys() and args[key]:
            runtime_parameters[key] = args[key]
            logging.debug('New value for {k}: {v}'.format(
                k=key, v=runtime_parameters[key]
            ))
    if args['primer3']:
        logging.debug('Found additional parameters for primer3.')
        for option, value in args['primer3']:
            try:
                value = ast.literal_eval(value)
                logging.debug('Evaluated {key} to {v}: {t}'.format(
                    key=option, v=value, t=type(value)
                ))
                primer3_options.update({option: value})
            except (SyntaxError, ValueError):
                logging.warning(
                    ('Could not parse option {key} with value {v} '
                     'from commandline.'
                     'Ignoring value').format(key=option, v=value)
                )
    logging.debug('Final parameter for primer3: {}'.format(primer3_options))



def main():
    """
    Delegates all tasks to the other functions.
    """
    # parse all arguments, let the argparse-module do its wonderful work
    args = parse_arguments()
    # setup logging
    setup_logging(args.loglevel)
    logging.debug('Received arguments: {}'.format(args))
    # is path to a config given?
    if args.config:
        runtime_parameters['config'] = args.config
    # create configparser object
    config = configparser.ConfigParser()
    # check whether path to configfile is a valid one
    if os.path.isfile(runtime_parameters['config']):
        # either default value or the one set per cmd argument is valid
        config.read(runtime_parameters['config'])
        logging.info(
            'Read config from file: {}'.format(runtime_parameters['config']))
    else:
        # no config read -> setting object to None
        config = None
        logging.info('No config passed to program.')
    parse_config_and_parameters(args, config)

    # # length of the insert between the primer pair
    # product_size_range = None
    # # are there any settings for primer3?
    # if config.has_section('primer3'):
    #     primer3_conf = config['primer3']
    #     logging.info('Found settings for primer3')
    #     logging.debug(
    #         'primer3-settings: {}'.format(' '.join(primer3_conf.keys())))
    #     if 'primer_product_size_range' in primer3_conf.keys():
    #         product_size_range = primer3_conf['primer_product_size_range']
    #
    # else:
    #     primer3_conf = None
    #     logging.warning(
    #         'Config contains no settings for primer3, '
    #         'running primer3 with default settings.')
    #
    # seq_included_region = extract_included_region(config['default'])
    #
    # """
    # Existing primer pairs specified via -p/--primerfiles will be read and bowtie
    # run against them.
    # """
    # if not args.keep_primer:
    #     logging.info('Generating new primer.')
    #     if args.additionalFasta is not None:
    #         sequences = args.additionalFasta
    #         logging.info(
    #             'Found additional file for primer generation {}'.format(
    #                 sequences.name
    #             ))
    #         logging.warning(
    #             'It is not possible to say whether a match found '
    #             'by bowtie is expected or not if an additional file for '
    #             'primer generation is passed.')
    #     else:
    #         sequences = args.FastaFile
    #         logging.info(
    #             'No additional file with sequences for primer-generation '
    #             'specified')
    #     """
    #     Extract sequence from FASTA file and get exact id of it as well
    #     """
    #     sequence, seq_id = parse_fasta(sequences, args.sequence)
    #     if not sequence:
    #         logging.error(
    #             "Could not find sequence with given ID-Prefix. Aborting")
    #         sys.exit(1)
    #     logging.info('Successfully extracted sequence')
    #
    #     """
    #     Generate primer pairs, depending on extracted sequence, primer3
    #     configuration and specified region for which primer shall be generated.
    #     """
    #     primer_dict = generate_primer(sequence, primer3_conf, args.primerfiles,
    #                                   seq_included_region)
    # else:
    #     logging.info('Trying to parse existing primer from files specified'
    #                  ' via -p/--primerfiles')
    #     logging.warning('Calculations whether a hit is expected or not '
    #                     'are based on the settings which would be used for '
    #                     'normal primer generation: SEQUENCE_INCLUDED_* for '
    #                     'the region, PRIMER_PRODUCT_SIZE_RANGE for size of '
    #                     'a possible insert and -s/--sequence for prefix-'
    #                     'matching with the reported id by bowtie in which '
    #                     'sequence the hit was found.')
    #     # store given sequence-parameter for prefix-matching with bowtie-results
    #     if args.sequence is not None:
    #         seq_id = args.sequence
    #     else:
    #         seq_id = None
    #
    #     primer_dict = parse_existing_primer(args.primerfiles)
    #     # we did not generate any primer, therefore we do not have any sequence
    #     # id to check whether a match is expected or not
    # # is an already existing bowtie-index specified?
    # if not args.index:
    #     # no index available, so we have to create our own one
    #     bowtie_index = setup_bowtie(args.FastaFile,
    #                                 args.loglevel == logging.DEBUG, args.bowtie)
    #     logging.info("No existing index for bowtie specified")
    # else:
    #     bowtie_index = args.index
    #     logging.info("Using existing bowtie-index")
    #
    # bowtie_result = run_bowtie(bowtie_index, args.primerfiles, args.bowtie,
    #                            args.loglevel == logging.WARNING,
    #                            product_size_range, args.bowtie_output)

    # create empty list for results
    # results = {}
    # for primer_tuple in bowtie_result:
    #     current_key, res = parse_bowtie_result(primer_tuple, config['default'],
    #                                            primer_dict, seq_included_region,
    #                                            args.additionalFasta is not None,
    #                                            seq_id, args.keep_primer)
    #     if res is not None:
    #         results.setdefault(current_key, []).append(res)
    #     else:
    #         continue
    #
    # # get value defined in config or fallback otherwise
    # max_number_matches = config['default'].getint('MAX_NUMBER_OF_MATCHES',
    #                                               MAX_NUMBER_OF_MATCHES)
    # logging.info('Maximal number of allowed matches set to {}'.format(
    #     max_number_matches))
    #
    # # set output either to STDOUT or output file if specified by user
    # if args.output is None:
    #     logging.info('No file for output specified, writing to STDOUT.')
    #     output = sys.stdout
    # else:
    #     logging.info('Output file specified, writing results to {}'.format(
    #         args.output.name))
    #     output = args.output
    #
    # # store intermediate all results which would be printed in output
    # printable_res = []
    # for key in sorted(results.keys()):
    #     matches = results[key]
    #     if len(matches) > max_number_matches:
    #         logging.debug(
    #             'Not printing results for {} because it has {} matches'.format(
    #                 key, len(matches)
    #             ))
    #     else:
    #         # add to intermediate results
    #         printable_res.append(matches)
    #
    # output.write(RESULT_HEADER + '\n')
    # for matches in sorted(printable_res, key=len):
    #     # write results
    #     output.write('\n'.join(matches))
    #     # final newline at end of results
    #     output.write('\n')


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

    def is_significant(values):
        # if last entry is a number it represents the number of matches at
        # the end of the primer; must not be lower than given value
        if values[-1].isdigit():
            if int(values[-1]) < bowtie_parse_options['LAST_MUST_MATCH']:
                logging.debug(
                    'Failed because last_must_match not fulfilled: {}'.format(
                        values))
                return False
            # only one number at the end and it is greate than last_must_match
            if len(values) == 1:
                return True
        # if last entry is a char there is no match, unless we do not care about
        # the entries at the end
        elif values[-1].isalpha() and bowtie_parse_options[
            'LAST_MUST_MATCH'] != 0:
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
        while bases_processed <= bowtie_parse_options['LAST_TO_CHECK']:
            current = values[-i]
            if current.isdigit():
                bases_processed += int(current)
            else:
                number_of_subs += 1
                bases_processed += 1
            i += 1
        if number_of_subs < bowtie_parse_options['LAST_MAX_ERROR']:
            return True
        else:
            logging.debug(
                'Failed because of {num} subs in last {last} '
                'bases, {max} allowed'.format(
                    num=number_of_subs,
                    last=bowtie_parse_options['LAST_TO_CHECK'],
                    max=bowtie_parse_options['LAST_MAX_ERROR']
                ))
        return number_of_subs < bowtie_parse_options['LAST_MAX_ERROR']

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


def default_string(key: str, dicts: dict) -> str:
    return " (default: " + str(dicts[key]) + ")"


def parse_arguments() -> argparse.Namespace:
    """
    This function parses the commandline arguments via the argparse-module,
    which additionally generates a help-hook, if a parameter is passed wrong.
    :return: A parser-object containing all parsed values.
    """
    parser = argparse.ArgumentParser(
        description="""
        genuprimer is able to generate new primer by calling primer3 and
        afterwards validate their uniqueness among other sequences by calling
        bowtie. The same can be accomplished for already existing primer pairs.
        """
    )
    parser.add_argument(
        "fasta_file", type=argparse.FileType('r'), metavar='path_to_fasta_file',
        help="File containing the sequences in FASTA-Format" + default_string(
            'fasta_file', runtime_parameters)
    )
    parser.add_argument(
        "-s", "--sequence", type=str, metavar='prefix_of_seq_id',
        help="""Partial ID of the sequence for which the primer
        shall be or have been generated.""", dest='seq_id'
    )
    parser.add_argument(
        '-c', '--config', type=str,
        help="""Configfile with various parameters. Has to
        include a '[default]'-section at top of the file, otherwise it can
        not be parsed.""" + default_string('config', runtime_parameters),
        metavar='path_to_config',
    )
    parser.add_argument(
        "-a", "--additionalFasta", type=argparse.FileType('r'),
        metavar='path_to_file', dest='additional_fasta',
        help="""An additional file containing the sequence for which primer
        shall be generated. First sequence inside of the file is taken if
        not specified otherwise via '-s'.""",
    )
    parser.add_argument(
        '--size', type=int, nargs=2, metavar=('min_size', 'max_size'),
        help="Size range of the product including primers."
    )
    parser.add_argument(
        '--pos', type=int, nargs=2, metavar=('begin', 'end'),
        help="Region between the primer which is not overlapped by them."
    )
    parser.add_argument(
        '-i', '--index', type=str,
        help="""If no bowtie-index is specified or found a new one will be generated
        for FastaFile. This option is directly forwarded to bowtie.""" +
             ' (default: ' + runtime_parameters['index'] + '/{FastaFile})',
    )
    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), nargs='?',
        help="""Output where the results should be stored.
        Default is standard output. Results are written as comma separated
        values (.csv).""" + ' (default: STDOUT)'
    )
    parser.add_argument(
        "--keep-primer", dest='keep_primer', action='store_true',
        help="""Set this option to start another run with the same primers from
        last run or some custom ones."""
    )
    parser.add_argument(
        '--last-must-match', dest='LAST_MUST_MATCH', type=int,
        help="How many of the last bases of a primer have to match to consider "
             "it a hit?" + default_string('LAST_MUST_MATCH',
                                          bowtie_parse_options)
    )
    parser.add_argument(
        '--last-to-check', dest='LAST_TO_CHECK', type=int,
        help="How many of the last bases of a primer should be checked "
             "considering LAST_MAX_ERROR." +
             default_string('LAST_TO_CHECK', bowtie_parse_options)
    )
    parser.add_argument(
        '--last-max-error', dest='LAST_MAX_ERROR', type=int,
        help="Maximum number of mismatches allowed to occur in the last "
             "LAST_TO_CHECK bases of a primer to consider it a hit." +
             default_string('LAST_MAX_ERROR', bowtie_parse_options)
    )
    parser.add_argument(
        '-l', '--limit-number-of-matches', dest='LIMIT_NUMBER_OF_MATCHES',
        help="Maximum number of hits of a primer pair before it is "
             "omitted from the results." +
             default_string(
                 'LIMIT_NUMBER_OF_MATCHES', bowtie_parse_options), type=int,
    )
    parser.add_argument(
        '--primer3', nargs=2, action='append',
        metavar=('primer3_option', 'value'), help="""Append any custom options
        for primer3 in a valid format for primer3-py. Options provided this way
        take precedence over values from configfile."""
    )
    parser.add_argument(
        "-p", "--primerfiles", type=str, help="""Prefix for the files where the
        primer pairs will be written to, if new ones are generated, or location
        of existing ones (see --keep-primer) with suffixes '_left.fas'
        and '_right.fas'.""" + default_string('prefix', runtime_parameters),
        metavar='prefix',
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
        "--show-bowtie", dest='show_bowtie_output', action='store_true',
        help="""Set this option to show the original results of bowtie, written
        to standard error output."""
    )
    parser.add_argument(
        "--bowtie", type=str, metavar='path_to_bowtie_executable',
        help="""The bowtie executable if not in PATH. If needed bowtie-build
        is expected to be found via appending '-build' to bowtie.""" +
             default_string('bowtie', runtime_parameters)
    )
    parser.set_defaults(keep_primer=False, show_bowtie_output=False)
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
