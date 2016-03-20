#!/usr/bin/env python3
import argparse

import primer3


def main():
    args = parse_arguments()
    sequence = parse_fasta(args.FastaFile, args.sequence)
    if not sequence:
        print("Could not find the sequence")
    find_primer(sequence)


def find_primer(sequence: str):
    primer3.setP3Globals({
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
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
    })
    print(res)


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
    return parser.parse_args()


def parse_fasta(fastaFile: argparse.FileType, seq_id: str):
    """
    Parses the submitted FASTA-File and extracts the sequence for primer3.
    :param seq_id:
    :param fastaFile:
    """
    seq = ""
    if not seq_id:
        fastaFile.readline()
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

    return seq


if __name__ == "__main__":
    main()
