# genuprimer

genuprimer is a tool to generate primer for a given sequence and afterwards looking for any matches of them against a desired
collection of sequences. It does so by using primer3 for generating the primer and bowtie for determining any matches.
This gives you the possibility to check to whether the primer pairs are only binding inside one sequence and at their
expected location only.

## Requirements

genuprimer uses and therefore needs [Primer3-py](https://libnano.github.io/primer3-py/index.html) to communicate with the 
primer3-library. A quick install guide can be found [here](https://libnano.github.io/primer3-py/quickstart.html#requirements).
This program has been developed for `Python3.4` but it may also work with newer versions.

## Options

A shortened and more compact version of available options and their meanings can also be found via `./genuprimer --help`.

Usage: 

`genuprimer.py [-h] [-s SEQUENCE] [-a ADDITIONALFASTA] [-p PRIMERFILES] [-i INDEX] [-d] [-v] -c CONFIG [-o OUTPUT]
                    [--bowtie BOWTIE] [--show-bowtie] [--keep-primer] FastaFile`
