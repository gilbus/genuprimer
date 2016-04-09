#!/usr/bin/env python3
import argparse
import os
import subprocess
from typing import Dict, Tuple

import primer3


def main():
    args = parse_arguments()
    sequence, sequence_id = parse_fasta(args.FastaFile, args.sequence)
    if not sequence:
        print("Could not find the sequence")
    primer_left, primer_right = find_primer(sequence)
    primerfile_left = open('primer_left.fas', 'w')
    primerfile_right = open('primer_right.fas', 'w')
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
    run_bowtie(args.FastaFile, bowtie_index)


def setup_bowtie(fasta_file: str) -> str:
    import re
    bowtie_index_dir = 'bowtie-index'
    os.makedirs(bowtie_index_dir, exist_ok=True)
    bowtie_index = "{index_dir}/{prefix}_bowtie".format(
        index_dir=bowtie_index_dir,
        prefix=
        re.split("/|\.", fasta_file.name)[-2])
    subprocess.run(
        "bowtie-build {fastaFile} {index_dir}".format(fastaFile=fasta_file.name,
                                                      index_dir=bowtie_index),
        shell=True)
    return bowtie_index


def run_bowtie(fasta_file: str, bowtie_index: str):
    subprocess.run(
        "bowtie -k 5000 -S -f {index} -1 primer_left.fas -2 primer_right.fas --sam-nohead".format(
            index=bowtie_index),
        shell=True)


def find_primer(sequence: str) -> Dict[str, str]:
    primer3.setP3Globals({
        'PRIMER_OPT_SIZE': 14,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 16,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_NUM_RETURN': 7,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[75, 100], [100, 125], [125, 150],
                                      [150, 175], [175, 200], [200, 225]],
    })

    sequence = sequence.replace('\n', '').replace('\r', '')
    res = primer3.bindings.designPrimers({
        'SEQUENCE_ID': 'mySequence',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': [36, 342]
    })  # type: dict
    primer_left = {}  # type: dict
    primer_right = {}  # type: dict
    for k in res.keys():
        line = k.split('_')
        if len(line) >= 4:
            if line[1] in ['RIGHT'] and \
                    line[3] == 'SEQUENCE':
                primer_right.update({k: res[k]})
            elif line[1] in ['LEFT'] and \
                    line[3] == 'SEQUENCE':
                primer_left.update({k: res[k]})
    return primer_left, primer_right


def parse_arguments():
    parser = \
        argparse.ArgumentParser(description=
                                """
                                Call primer3 to generate primers and validate
                                their uniqueness among other
                                sequences with the help of bowtie
                                """
                                )
    parser.add_argument("FastaFile", type=argparse.FileType('r'),
                        help="File containing the sequences in FASTA-Format")
    parser.add_argument("-s", "--sequence", type=str,
                        default=None,
                        help=
                        """
                        ID of the sequence which should be given to primer3
                        (default: first sequence in FASTA-File).
                        If omitted prefix-matching is used for identification
                        and first hit will be used.
                        """
                        )
    parser.add_argument('-i', '--index', type=str, default=None,
                        help="""Use existing bowtie-index. This option is
                        directly forwarded to bowtie."""
                        )
    return parser.parse_args()


def parse_fasta(fastaFile: argparse.FileType, seq_id: str) -> Tuple[str, str]:
    """
    Parses the submitted FASTA-File and extracts the sequence for primer3.
    :param seq_id:
    :param fastaFile:
    """
    seq = ""
    if not seq_id:
        seq_id = fastaFile.readline()[1:]
        for line in fastaFile:
            if line in ['\n', '\r\n']:
                break
            seq += line
    else:
        for line in fastaFile:
            if line[0] == '>' and line.startswith(seq_id, 1):
                break
        for line in fastaFile:
            if line in ['\n', '\r\n']:
                break
            seq += line

    return seq, seq_id


if __name__ == "__main__":
    main()
