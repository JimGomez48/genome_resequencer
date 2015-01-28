
NAME
=====
Genome Resequencer


AUTHORS
========
James Gomez


DESCRIPTION
============
This software is a set of python scripts for resequencing a donor genome and
finding SNP and indel variants with respect to a given reference genome. It
implements a trivial exhaustive search alignment algorithm as a baseline
approach and a hashtable-based approach where a k'mer index of the reference
genome is built to search for partial matches of reads in constant time. Example
genomes are included.


REQUIREMENTS
=============
- Python 2.7.x
- NumPy numerical and scientific computing package for Python


USAGE
======
python <name_of_script> <arguments> <options>

The following scripts are included and should be run in the specified order to
perform proper resequencing:

1) align_reads.py:
        Aligns a set of reads from a donor genome to a reference genome. Outputs
        the aligned reads as a separate file 'align_*.txt'. Use -h option for
        more details.

1) pile_up.py:
        Generates a pile-up from the aligned reads and then generates a
        consensus donor genome from the pile-up. Outputs the consensus donor as
        a separate file 'cons_*.txt'. Use -h option for more details.

3) find_variants.py:
        Finds SNPs and indels in the donor genome with respect to the reference
        and outputs the result as a file 'variants_*.txt'. Use -h option for
        more details.

* For more information on how to run a script, use the -h option

Also included are the following support scripts and files:

- edit_distance.py:
        A support file containing methods which implement various edit-distance
        algorithms. This file's methods should not be called directly.

- common.py:
        A support file containing convenience methods. This file's methods
        should not be called directly.

- genomes.zip:
        A set of genomes, 10k, 100k, 1m, 10m, each of which has the following
        files:
            - ref_*.txt: contains the reference genome sequence
            - donor_*.txt: contains the donor genome sequence
            - reads_*.txt: contains the set of reads from the donor genome
            - ans_*.txt: contains the variants in donor with respect to its reference
        The sample donor genomes contain only SNP and indel variants. Reads have
        an individual base error rate of 0.01 and garbage read rate of 0.05.