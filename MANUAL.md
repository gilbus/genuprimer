# genuprimer

genuprimer is a tool to generate primer for a given sequence and afterwards looking for any matches
of them against a desired collection of sequences. It does so by using primer3 for generating the
primer and bowtie for determining any matches. This gives you the possibility to check to whether
the primer pairs are only binding inside one sequence and at their expected location only.

## Requirements

genuprimer uses and therefore needs [Primer3-py](https://libnano.github.io/primer3-py/index.html)
to communicate with the primer3-library. A quick install guide can be found 
[here](https://libnano.github.io/primer3-py/quickstart.html#requirements).
This program has been developed for `Python3.4` but it may also work with newer versions.
Since primer3 is only needed for primer generation, which can be skipped via `--keep-primer`, it
can still be used for primer validation. Take a look at the Examples-Section below.

## Options

A shortened and more compact version of available options and their meanings can also be 
found via `./genuprimer --help`.

Usage: 

`genuprimer.py [-h] [-s prefix_of_seq_id] [-c path_to_config]
                     [-a path_to_file] [--size min_size max_size]
                     [--pos begin end] [-i INDEX] [-o [OUTPUT]]
                     [--keep-primer] [--last-must-match LAST_MUST_MATCH]
                     [--last-to-check LAST_TO_CHECK]
                     [--last-max-error LAST_MAX_ERROR]
                     [-l LIMIT_NUMBER_OF_MATCHES]
                     [--primer3 primer3_option value] [-p prefix] [-v] [-d]
                     [--show-bowtie] [--bowtie path_to_bowtie_executable]
                     path_to_fasta_file`

### Positional arguments

  `path_to_fasta_file`
        File containing the sequences in valid FASTA format. In the rest of the manual it will be
        referred to as `FastaFile`.
 

### Optional arguments

  `-h, --help`            show this help message and exit

  `-s prefix_of_seq_id, --sequence prefix_of_seq_id`
        Partial ID of the sequence for which the primer shall be or have been generated. To identify
        the correct FASTA sequence prefix matching is applied and the first sequence where a match
        is found is read from the file. If no primer generation is performed, see `--keep-primer`,
        the supplied value is used for bowtie to see whether a possible match is expected or not.
        If no value is supplied the **first sequence found** is taken.

  `-c path_to_config, --config path_to_config`
        Configfile with various parameters concerning the primer generation and the evaluation of
        the bowtie result. See the Config-Section below for examples and instructions concerning the
        format.
        If a file `genuprimer.conf` is found and readable in the same directory it will be used.

  `-a path_to_file, --additionalFasta path_to_file`
        An additional FASTA file containing the sequence for which primer shall be generated. If
        a file is omitted `-s prefix_of_seq_id` is used to identify the chosen sequence.
        Additionally it is expected that the sequence chosen from `additionalFasta` is **not**
        present in `FastaFile`. Since `bowtie` is applied against `FastaFile` every hit will
        be marked as not expected.

  `--size min_size max_size`
        Size range of the product including primers. See Examples-Section below.

  `--pos begin end`
        Positions of the region of interest between the left and the right primer which is not 
        overlapped by them. Must be absolute Positions inside the chosen Sequence.

  `-i INDEX, --index INDEX`
        If no bowtie-index is specified or found at the default location a new one will be
        generated for FastaFile. This option is directly forwarded to bowtie. 
        The default location is `bowtie-index/{generated_name}` where `generated_name` is calculated
        by taking the value of `FastaFile`, splitting it at every occurrence of `/` or `.` and the
        second last value is taken and `_bowtie` appended to it. See Examples-Section below.

  `-o [OUTPUT], --output [OUTPUT]`
        Output where the final results shall be stored. Default is STDOUT, e.g. printing to the
        console. If an existing file is specified its contents will be overwritten.
        Results are written as comma separated values (.csv). See Result-Section below.

  -p prefix, --primerfiles prefix
                        Prefix for the files where the primer pairs will be
                        written to, if new ones are generated, or location of
                        existing ones (see --keep-primer) with suffixes
                        '_left.fas' and '_right.fas'. (default: genuprimer)

  `--keep-primer`
        If this option is omitted no new primers will be generated. Instead, existing primer pairs,
        again located via `prefix` in combination with the suffixes, will be read and used for
        `bowtie`. This option can be appended to any previous command where primer have been
        generated but maybe you want to change the following settings concerning the evaluation of 
        the `bowtie` results to get different final results.
        Useful if you want to check some existing primer pairs or the generation takes very long.
        If one of the files contains more primers than the other one the 'lonely' ones at the end of
        the larger file will be dropped.

  `--last-must-match LAST_MUST_MATCH`
        How many of the last bases of a primer have to match to consider it a hit? (default: 3)

  --last-to-check LAST_TO_CHECK
                        How many of the last bases of a primer should be
                        checked considering LAST_MAX_ERROR. (default: 12)
  --last-max-error LAST_MAX_ERROR
                        Maximum number of mismatches allowed to occur in the
                        last LAST_TO_CHECK bases of a primer to consider it a
                        hit. (default: 5)
  -l LIMIT_NUMBER_OF_MATCHES, --limit-number-of-matches LIMIT_NUMBER_OF_MATCHES
                        Maximum number of hits of a primer pair before it is
                        omitted from the results. (default: 5)
  --primer3 primer3_option value
                        Append any custom options for primer3 in a valid
                        format for primer3-py. Options provided this way take
                        precedence over values from configfile.
  -v, --verbose         Be verbose by showing INFO messages.
  -d, --debug           Print lots of DEBUG messages.
  --show-bowtie         Set this option to show the original results of
                        bowtie, written to standard error output.
  --bowtie path_to_bowtie_executable
                        The bowtie executable if not in PATH. If needed
                        bowtie-build is expected to be found via appending
                        '-build' to bowtie. (default: bowtie)
