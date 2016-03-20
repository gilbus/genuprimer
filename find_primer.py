#!/usr/bin/env python3
import argparse


def main():
    args = parse_arguments()
    sequence = parse_fasta(args.FastaFile, args.sequence)
    print(sequence)


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


def parse_fasta(fastaFile: argparse.FileType, sequence: str):
    """
    Parses the submitted FASTA-File and extracts the sequence for primer3.
    """
    seq = ""
    if not sequence:
        fastaFile.readline()
        for line in fastaFile:
            if line in ['\n', '\r\n']:
                break
            seq += line
    else:
        for line in fastaFile:
            if line
    return seq



if __name__ == "__main__":
    main()
