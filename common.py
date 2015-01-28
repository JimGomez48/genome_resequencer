from sys import stdout
from math import ceil


REF_PRE = 'ref_'
READS_PRE = 'reads_'
DONOR_PRE = 'donor_'
ANS_PRE = 'ans_'
MYANS_PRE = 'variants_'
ALIGN_PRE = 'align_'
CONS_PRE = 'cons_'

# REFS_DIR = 'refs/'
# READS_DIR = 'reads/'
# DONORS_DIR = 'donors/'
# ANS_DIR = 'ans/'
# MYANS_DIR = 'myans/'
# ALIGN_DIR = 'alignments/'
# CONS_DIR = 'cons/'

REFS_DIR = READS_DIR = DONORS_DIR = ANS_DIR = MYANS_DIR = ALIGN_DIR = CONS_DIR = ''


def print_progress(progress, msg='Progress'):
    """
    Decimal progress out of 1.0
    :param progress: the progress value [0.0, 1.0]
    :param msg: the progress message to prepend to the progress output
    """
    progress = min(100., progress * 100.)
    stdout.write("\t%s... [%.2f%%]\r" % (msg, progress))
    stdout.flush()


def load_genome(filename):
    """
    Loads a genome from file and returns the genome as a single string
    :param filename: the file containing the genome
    :return: the genome string and the name of the genome as a tuple
    """
    genome = ''
    with open(filename) as ref_file:
        line = ref_file.next()
        if not '>' in line:
            raise Exception('Could not find genome name')
        else:
            name = str(line.rstrip())
            name = name[1: len(name)]
            # print "\tLoading genome from '" + filename + "'"
            for line in ref_file:
                if '>' in line:
                    continue
                line = line.strip()
                genome += line
    return genome, name


def compute_kmer_length(read_length, max_mismatches):
    """
    Compute the optimal k'mer length for the given read length and desired
    maximum mismatch
    :param read_length: the average length of the read
    :param max_mismatches: the max number of mismatches allowed during alignment
    :return: the optimal k'mer length
    """
    num_segments = max_mismatches + 1
    while 1:
        if read_length % num_segments == 0:
            break
        num_segments += 1
    return int(ceil(float(read_length) / float(num_segments)))