import argparse
import numpy as np
import common as cm

def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="Aligns reads from a donor genome to a given reference "
                    "genome. Produces an alignment file as output"
    )
    arg_parser.add_argument(
         'ref_file',
         type=str,
         help='the file containing the reference genome')
    arg_parser.add_argument(
         'reads_file',
         type=str,
         help='the file containing reads from the donor genome')
    # arg_parser.add_argument(
    #     '-g', '--genome-name',
    #     metavar='CHOICE',
    #     required=False,
    #     type=int,
    #     choices=range(1, 4),
    #     default=1,
    #     help='(1) hw2_grad (2) practice_E_1 (3) practice_W_3. Default: 1')
    # arg_parser.add_argument(
    #     '-m', '--max-mismatch',
    #     metavar='NUM',
    #     required=False,
    #     type=int,
    #     default=2,
    #     help='max number of mismatches allowed per read. Default: 2')
    arg_parser.add_argument(
        '-t', '--trivial-align',
        required=False,
        action='store_true',
        default=False,
        help='use the trivial algorithm to align reads to the reference')
    args = arg_parser.parse_args()
    print_args(args)
    return args


def print_args(args):
    print '=========================================='
    print 'ref-file:\t' + str(args.ref_file)
    print 'reads-file:\t' + str(args.reads_file)
    if args.trivial_align:
        print 'using index:\tFALSE'
    else:
        print 'using index:\tTRUE'
    # print 'max-mismatch:\t' + str(args.max_mismatch)
    print '=========================================='


def create_genome_index(genome, k=10):
    print 'Creating ref-genome index...'
    index = {}
    for i in xrange(0, len(genome) - k + 1):
        sequence = ''.join(genome[i: i + k])
        index.setdefault(sequence, []).append(int(i))
    print '\tsize: ' + str(len(index)) + ' keys'
    return index


def align_read_trivial(ref, read, k):
    best_pos = None
    min_mismatch = 100
    thresh = 2
    # for every starting position in the reference genome
    for i in xrange(len(ref) - len(read) + 1):
        try:
            num_mismatches = sum(1 for j in xrange(len(read)) if ref[i+j] != read[j])
        except IndexError as e:
            print i
            print len(ref)
            raise e
        if num_mismatches < min_mismatch:
            min_mismatch = num_mismatches
            best_pos = i
            if num_mismatches <= thresh:
                break
    return best_pos, min_mismatch


def align_read_hashtable(ref, read, index, k):
    best_pos = None
    min_mismatch = 100
    # try to hash all distinct k'mers of the read
    for start in xrange(0, len(read), k):
        if start + k > len(read):            start = len(read) - k
        kmer = read[start: start + k]
        try:
            positions = index[kmer]
            for p in positions:
                try:
                    # get number of mismatches between the read and reference
                    # at the hashed position
                    n_mismatches = 0
                    for j in xrange(len(read)):
                        if read[j] != ref[p + j - start]:
                            n_mismatches += 1
                except IndexError:
                    continue
                # if this position results in fewer mismatches, save the new
                # position and min-mismatch values
                if n_mismatches < min_mismatch:
                    min_mismatch = n_mismatches
                    best_pos = p - start
        except KeyError:  # couldn't hash -> k'mer not in the reference
            continue
    return best_pos, min_mismatch


def map_reads_forward(ref, ref_name, reads_file, k, trivial=False):
    if not trivial:
        index = create_genome_index(ref, k)
    align_file = cm.ALIGN_DIR + cm.ALIGN_PRE + ref_name + '.txt'
    print 'Mapping reads...'
    with open(reads_file, 'r') as reads:
        num_lines = float(sum(1 for line in reads))
    with open(reads_file, 'r') as reads, open(align_file, 'w') as outfile:
        outfile.write('>' + ref_name + '\n')
        line_num = 0.
        for line in reads:
            if '>' in line:
                continue
            # get reads from the line (already oriented correctly)
            read1, read2 = line.strip().split(',')
            # map the reads
            if trivial:
                best_pos_1, min_mismatch_1 = align_read_trivial(ref, read1, k)
                best_pos_2, min_mismatch_2 = align_read_trivial(ref, read2, k)
            else:
                best_pos_1, min_mismatch_1 = align_read_hashtable(ref, read1, index, k)
                best_pos_2, min_mismatch_2 = align_read_hashtable(ref, read2, index, k)
            # write the aligned reads to file
            outfile.write(read1 + ':' + str(best_pos_1) + ',' +
                          read2 + ':' + str(best_pos_2) + '\n')
            line_num += 1
            # cm.print_progress(lin_num / 30000000.0)
            cm.print_progress(line_num / num_lines)
    print '\tComplete'
    return align_file


def map_reads_inverted(ref, ref_name, reads_file, k):
    index = create_genome_index(ref, k)
    align_file = cm.ALIGN_DIR + cm.ALIGN_PRE + ref_name + '.txt'
    print 'Mapping reads...'
    with open(reads_file, 'r') as reads:
        num_lines = float(sum(1 for line in reads))
    print num_lines
    with open(reads_file, 'r') as reads, open(align_file, 'w') as outfile:
        outfile.write('>' + ref_name + '\n')
        line_num = 0.
        for line in reads:
            if '>' in line:
                continue
            # get reads from the line and orient them in the same direction
            read1, read2 = line.strip().split(',')
            read2 = read2[::-1]
            # map the reads in the current orientation
            best_pos1_f, min_mismatch1_f = align_read_hashtable(ref, read1, index, k)
            best_pos2_f, min_mismatch2_f = align_read_hashtable(ref, read2, index, k)
            # reverse the reads
            read1 = read1[::-1]
            read2 = read2[::-1]
            # map the reads in the reversed orientation
            best_pos1_r, min_mismatch1_r = align_read_hashtable(ref, read1, index, k)
            best_pos2_r, min_mismatch2_r = align_read_hashtable(ref, read2, index, k)
            # if the better mapping was in the first orientation, flip the
            # reads back to the original orientation and write their alignment
            if min_mismatch1_f < min_mismatch1_r:
                # flip the reads back to original orientation
                read1 = read1[::-1]
                read2 = read2[::-1]
                outfile.write(read1 + ':' + str(best_pos1_f) + ',' +
                          read2 + ':' + str(best_pos2_f) + '\n')
            else:
                outfile.write(read1 + ':' + str(best_pos1_r) + ',' +
                          read2 + ':' + str(best_pos2_r) + '\n')
        line_num += 1
        # cm.print_progress(lin_num / 30000000.0)
        cm.print_progress(line_num / num_lines)
    print '\tComplete'
    return align_file


def pretty_print_aligned_reads(ref, align_file):
    reads = []
    alignments = []
    # create list of spaced alignments and list of paired-read position tuples
    with open(align_file, 'r') as infile:
        for line in infile:
            alignment1, alignment2 = line.strip().split(',')
            read1, position1 = alignment1.split(':')
            read2, position2 = alignment2.split(':')
            if position1 == 'None' or position2 == 'None':
                continue
            space = '.' * (int(position2) - int(position1) - len(read1))
            reads.append(read1 + space + read2)
            alignments.append((int(position1), int(position2)))
    # remove bad alignments from the reads
    good_alignments = [120 < x[1] - x[0] < 180 for x in alignments]
    best_reads = [reads[i] for i in range(len(good_alignments))
                  if good_alignments[i]]
    best_alignments = [alignments[i] for i in range(len(alignments))
                       if good_alignments[i]]
    del good_alignments
    # create list of sorted indeces to alignments based on first read position
    first_alignment = [x[0] for x in best_alignments]
    alignment_indices = np.argsort(first_alignment)
    sorted_reads = [best_reads[i] for i in alignment_indices]
    sorted_alignments = [best_alignments[i] for i in alignment_indices]
    # print the alignments against the reference genome
    active_reads = []
    line_length = 100
    print '\n\n' + '=' * (line_length + 6) + '\n'
    for i in range(len(ref) / line_length):
        next_ref = ref[i * line_length: (i + 1) * line_length]
        new_read_indices = [j for j in range(len(sorted_reads))
                            if i * line_length <=
                            sorted_alignments[j][0] <
                            (i + 1) * line_length]
        space_amounts = [sorted_alignments[index][0] %
                         line_length for index in new_read_indices]
        new_reads = [sorted_reads[index] for index in new_read_indices]
        new_reads_spaces = [' ' * space_amounts[j] + new_reads[j]
                                 for j in range(len(new_reads))]
        empty_active_read_indices = [index for index in range(len(active_reads))
                                     if active_reads[index] == '']
        for j in range(min(
                len(new_reads_spaces), len(empty_active_read_indices))):
            active_reads[empty_active_read_indices[j]] = new_reads_spaces[j]

        if len(new_reads_spaces) > len(empty_active_read_indices):
            active_reads += new_reads_spaces[len(empty_active_read_indices):]
        printed_reads = ['Read: ' + read[:line_length] for read in active_reads]
        active_reads = [read[line_length:] for read in active_reads]
        while len(active_reads) > 0:
            last_thing = active_reads.pop()
            if last_thing != '':
                active_reads.append(last_thing)
                break
        output_lines = ['Ref:  ' + next_ref] + printed_reads
        output_str = 'Reference index: ' + str(i * line_length) + \
                     '\n' + '\n'.join(output_lines) + '\n\n' + \
                     '-' * (line_length + 8) + '\n'
        print output_str


def main():
    args = parse_args()

    # get the ref and reads file paths
    ref_file = args.ref_file
    reads_file = args.reads_file

    # load the reference genome as a string
    ref, ref_name = cm.load_genome(ref_file)
    # kmer_length = cm.compute_kmer_length(50, args.max_mismatch)

    # align the reads to the reference
    # map_reads_inverted(ref, ref_name, reads_file, k=10, trivial=args.trivial_align)
    map_reads_forward(ref, ref_name, reads_file, k=10, trivial=args.trivial_align)

    # pretty_print_aligned_reads(ref, align_file)
    print 'DONE'


if __name__ == '__main__':
    main()