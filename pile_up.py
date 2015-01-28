import argparse
import numpy as np
import common as cm


def __parse_args():
    arg_parser = argparse.ArgumentParser(
        description="Generates a consensus donor genome with respect to "
                    "a reference genome from aligned reads."
    )
    arg_parser.add_argument(
         'ref_file',
         type=str,
         help='the file containing the reference genome')
    arg_parser.add_argument(
         'align_file',
         type=str,
         help='the file containing the aligned reads from the donor genome')
    args = arg_parser.parse_args()
    __print_args(args)
    return args


def __print_args(args):
    print '=========================================='
    # print 'genome-name:\t' + str(cm.get_genome_name(args.genome_name))
    print 'ref-file:\t' + str(args.ref_file)
    print 'align-file:\t' + str(args.align_file)
    print '=========================================='


def __get_consensus(ref, ref_name, align_file):
    print 'Generating pile-up...'
    index_allele_map = ['A', 'C', 'G', 'T']
    allele_index_map = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
    }
    NONE = 'None'
    # get consensus pile_up
    consensus_pile = np.zeros((len(ref), 4), dtype=np.uint32)
    with open(align_file, 'r') as alignment_file:
        num_lines = float(sum((1 for line in alignment_file if not '>' in line)))
    with open(align_file, 'r') as alignment_file:
        line_count = 1.0
        for line in alignment_file:
            if '>' in line:
                continue
            alignment1, alignment2 = line.strip().split(',')
            read1, position1 = alignment1.split(':')
            read2, position2 = alignment2.split(':')
            # if paired reads weren't mapped, continue
            if position1 == NONE or position2 == NONE:
                continue
            # convert positions to ints
            position1 = int(position1)
            position2 = int(position2)
            # if paired reads weren't mapped, continue
            if position1 == -1 or position2 == -1:
                continue
            # if spacing between reads out of whack
            # don't contribute to pile-up
            MAX_GAP = 50
            if position2 - position1 < 100 or position2 - position1 > 200:
                continue
            # contribute read1 to consensus pile
            for i in xrange(len(read1)):
                c = read1[i]
                consensus_pile[position1 + i][allele_index_map[c]] += 1
            # contribute read2 to consensus pile
            for i in xrange(len(read2)):
                c = read2[i]
                consensus_pile[position2 + i][allele_index_map[c]] += 1
            line_count += 1.0
            progress = line_count / num_lines
            cm.print_progress(progress)
    cm.print_progress(1)
    print '\tComplete'
    # get consensus string from pileup
    print 'Resolving consensus donor...'
    consensus_file_name = cm.CONS_DIR + cm.CONS_PRE + ref_name + '.txt'
    with open(consensus_file_name, 'w') as consensus_file:
        consensus_file.write('>' + ref_name + '\n')
        for i in xrange(len(consensus_pile)):
            if i > 0 and i % 80 == 0:
                consensus_file.write('\n')
            msum = sum(consensus_pile[i])
            max_pos = np.argmax(consensus_pile[i])
            mmax = consensus_pile[i][max_pos]
            # if it's a true max, write the consensus allele
            if msum > 0 and float(mmax) / float(msum) > 0.25:
                consensus_file.write(index_allele_map[max_pos])
            else:  # otherwise, stick with the ref allele (too much ambiguity)
                consensus_file.write(ref[i])
            cm.print_progress(float(i) / float(len(consensus_pile)))
    cm.print_progress(1)
    print '\tComplete'
    return consensus_file_name


def pile_up(ref_genome, ref_name, align_file_name):
    return __get_consensus(ref_genome, ref_name, align_file_name)


def pretty_print_ref_consensus(ref, consensus_file):
    # load consensus string
    consensus = ''
    with open(consensus_file, 'r') as c_file:
        for line in c_file:
            if '>' in line:
                continue
            consensus += line.strip()
    # print consensus against ref
    line_length = 100
    for i in xrange(0, len(ref), line_length):
        print '-' * (line_length + 6)
        print
        print 'ref-pos: ' + str(i)
        print 'ref:  ' + ref[i: i + min(line_length, len(ref) - i)]
        match_line = '      '
        for j in xrange(min(line_length, len(consensus) - i)):
            if ref[i + j] != consensus[i + j]:
                match_line += '*'
            else:
                match_line += ' '
        print match_line
        print 'cons: ' + consensus[i: i + min(line_length, len(consensus) - i)]
        print


def __main():
    args = __parse_args()

    # get the ref and consensus file paths
    ref_file = args.ref_file
    align_file = args.align_file

    ref, ref_name= cm.load_genome(ref_file)

    consensus_file = pile_up(ref, ref_name, align_file)
    # pretty_print_ref_consensus(ref, consensus_file)
    print 'DONE'


if __name__ == '__main__':
    __main()
