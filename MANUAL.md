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

`genuprimer.py [-h] [-s prefix_of_seq_id] [-c path_to_config][-a path_to_file]
[--size min_size max_size][--pos begin end] [-i INDEX] [-o [OUTPUT]][--keep-primer]
[--last-must-match LAST_MUST_MATCH][--last-to-check LAST_TO_CHECK]
[--last-max-error LAST_MAX_ERROR][-l LIMIT_NUMBER_OF_MATCHES][--primer3 primer3_option value] 
[-p prefix][-v][-d][--show-bowtie] [--bowtie path_to_bowtie_executable]
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

  `-p prefix, --primerfiles prefix`
        Prefix for the files where the primer pairs will be written to, if new ones are generated,
        or location of existing ones, see next option, with suffixes '_left.fas' and '_right.fas'.
        See Examples-Section below.

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
        How many of the last bases of a primer have to match to consider it a hit? See
        Bowtie-Section below. The default value is 3.

  `--last-to-check LAST_TO_CHECK`
        How many of the last bases of a primer should be checked considering LAST_MAX_ERROR, see 
        next option.See Bowtie-Section below. The default value is 12.

  `--last-max-error LAST_MAX_ERROR`
        Maximum number of mismatches allowed to occur in the last LAST_TO_CHECK bases of a primer to
        consider it a hit.See Bowtie-Section below. The default value is 5.

  `-l LIMIT_NUMBER_OF_MATCHES, --limit-number-of-matches LIMIT_NUMBER_OF_MATCHES`
        Maximum number of hits of a primer pair before it is not shown in the results. The default
        value is 5.

  `--primer3 primer3_option value`
        Append any custom options for primer3 in a valid format for primer3-py. Options provided 
        this way take precedence over values from configfile. The values must be correct formatted
        in python notation since they are directly forwarded to primer3-py.

  `-v, --verbose`
        Be verbose by showing INFO messages. These are just indicators what the program is
        currently doing and activating them is highly recommended.

  `-d, --debug`
        Print lots of DEBUG messages. Useful if you do not understand actions of the program and why
        it is executing them or how your config or command line arguments are parsed and evaluated.

  `--show-bowtie`
        Set this option to show the original results/output of bowtie, written to standard error 
        output (STDERR).

  `--bowtie path_to_bowtie_executable`
        The bowtie executable if not in PATH. If needed bowtie-build is expected to be found via 
        appending '-build' to the bowtie-call.

## Config
A config file can be passed via `-c path_to_config`, see the beginning of the upper section.
It is parsed by the [configparser](https://docs.python.org/3.4/library/configparser.html) and
follows INI
[syntax](https://docs.python.org/3.4/library/configparser.html#supported-ini-file-structure).
It consists of two possible sections.

### default
Available options are:
* `TARGET_POSITION_BEGIN`
    Position of the begin of the region of interest, see `--pos`
* `TARGET_POSITION_END`
    End of the begin of the region of interest, see `--pos`
* `PRIMER_PRODUCT_SIZE_MIN`
    Minimal size of the product including the sizes of the primer, see `--size`
* `PRIMER_PRODUCT_SIZE_MAX`
    Maximal size of the product including the sizes of the primer, see `--size`
* `LAST_MUST_MATCH`
    Shortcut for same named command line option `--last-must-match`
* `LAST_TO_CHECK`
    Shortcut for same named command line option `--last-to-check`
* `LAST_MAX_ERROR`
    Shortcut for same named command line option `--last-max-error`

### primer3
Options which are passed to primer3-py containing criteria for the primer generation.

#### Example
The following example config will also be used in the Examples-Section

    [default]
    # These four values are mandatory and must be specified either via command line or here
    TARGET_POSITION_BEGIN = 4709738
    TARGET_POSITION_END   = 4710138
    PRIMER_PRODUCT_SIZE_MIN = 500
    PRIMER_PRODUCT_SIZE_MAX = 700

    # These three values are having default values
    LAST_MUST_MATCH = 3
    LAST_TO_CHECK = 12
    LAST_MAX_ERROR = 5
    
    # This section must not be defined necessarily but it is the only way to influence primer
    # generation
    [primer3]
    PRIMER_OPT_SIZE = 20
    PRIMER_INTERNAL_MAX_SELF_END = 8
    PRIMER_MIN_SIZE = 10
    PRIMER_MAX_SIZE = 30
    PRIMER_OPT_TM = 60
    PRIMER_MIN_TM = 57.0
    PRIMER_MAX_TM = 63.0
    PRIMER_MIN_GC = 20.0
    PRIMER_MAX_GC = 80.0
    PRIMER_MAX_POLY_X = 100
    PRIMER_NUM_RETURN = 10
    PRIMER_INTERNAL_MAX_POLY_X = 100
    PRIMER_SALT_MONOVALENT = 50.0
    PRIMER_DNA_CONC = 50.0
    PRIMER_MAX_NS_ACCEPTED = 1
    PRIMER_MAX_SELF_ANY = 12
    PRIMER_MAX_SELF_END = 8
    PRIMER_PAIR_MAX_COMPL_ANY = 12
    PRIMER_PAIR_MAX_COMPL_END = 8
    PRIMER_LIBERAL_BASE = 1

Again, using a config file is **not** mandatory but simplifies the whole process.

## Result
By using the config above with `PRIMER_NUM_RETURN = 10` ten primer pairs will be generated and
possible matches found by bowtie are evaluated. **Not all matches found by bowtie are included
here**, since these are evaluated, see the Bowtie-Section for more information. If a match is
detected to be considerable it can be found here and all matches are also marked *expected* or 
*not expected*. A match is considered expected if its position, size and name of the sequence where
the match has been found comply with the settings used to generate the primer pairs.

The final results will be written as comma-separated-values (csv).

    FWD_ID,REV_ID,MATCH_ID,FWD,REV,START,STOP,LENGTH,EXP
    PRIMER_LEFT_0_SEQUENCE,PRIMER_RIGHT_0_SEQUENCE,Chr1,GTGGTATTGCGTTCGCTTCG,TGGTGACTTAAGCGACTTGC,4709649,4710333,684,1

One considerable match is represented by nine values:
* `FWD_ID`
    ID of the forward primer where the match has been found
* `REV_ID`
    bowtie is configured to look for matches where the left and right primer are in correct distance
    to each other, see next section for further details. Therefore, one reverse primer also has to 
    match..
* `MATCH_ID`
    ID of the sequence inside `FastaFile` where the match has been found.
* `FWD`
    Sequence of the forward match.
* `REV`
    Sequence of the reverse match.
* `START`
    Begin position of the match inside `FastaFile`
* `END`
    End position of the match inside `FastaFile`.
* `LENGTH`
    Length of the whole match flanked by the primer pair.
* `EXP`
    Whether a match is considered expected, see above explanation.

## Bowtie
### Call
bowtie is configured to consider the left and right primer as 'paired-end reads' and the `--size`
values are also forwarded as `-I/--minins` and `-X/--maxins`.

### Evaluation
## Examples
