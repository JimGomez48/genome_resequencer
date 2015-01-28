import argparse
import common as cm
import edit_distance as ed


def __parse_args():
    arg_parser = argparse.ArgumentParser(
        description="Finds variants in a consensus donor genome with respect to "
                    "a reference genome."
    )
    arg_parser.add_argument(
         'ref_file',
         type=str,
         help='the file containing the reference genome')
    arg_parser.add_argument(
         'cons_file',
         type=str,
         help='the file containing the consensus donor genome')
    # arg_parser.add_argument(
    #     '-g', '--genome-name',
    #     metavar='CHOICE',
    #     required=False,
    #     type=int,
    #     choices=range(1, 4),
    #     default=1,
    #     help='(1) hw2_grad (2) practice_E_1 (3) practice_W_3. Default: 1')
    arg_parser.add_argument(
        '-e', '--edit-dist',
        metavar='STRATEGY',
        required=False,
        type=str,
        choices=['a', 'g', 'l'],
        default='a',
        help='(a) affine-gap alignment (g) global alignment (l) local aligment. Default: \'a\'')
    # arg_parser.add_argument(
    #     '-t', '--snp-thresh',
    #     metavar='T',
    #     required=False,
    #     type=int,
    #     default=2,
    #     help='max number of mismatches allowed in window')
    args = arg_parser.parse_args()
    __print_args(args)
    return args


def __print_args(args):
    print '=========================================='
    # print 'genome-name:\t' + str(cm.get_genome_name(args.genome_name))
    print 'ref-file:\t' + str(args.ref_file)
    print 'cons-file:\t' + str(args.cons_file)
    # print 'snp-thresh:\t' + str(args.snp_thresh)
    if args.edit_dist == 'a':
        print 'alignment:\taffine-gap'
    elif args.edit_dist == 'g':
        print 'alignment:\tglobal'
    elif args.edit_dist == 'l':
        print 'alignment:\local'
    print '=========================================='


def __write_answer_file(ref_name, snps, inserts, deletes, edist_strategy):
    print 'Writing answer file...'
    # sort variants by position in the reference
    print '\tsorting answers...'
    if len(snps) > 0:
        snps.sort(key=lambda x: x[2])
    if len(inserts) > 0:
        inserts.sort(key=lambda x: x[1])
    if len(deletes) > 0:
        deletes.sort(key=lambda x: x[1])
    # write variants to file
    ans_file = cm.MYANS_DIR + cm.MYANS_PRE
    if edist_strategy == 'a':
        ans_file += 'affine_' + ref_name + '.txt'
    elif edist_strategy == 'g':
        ans_file += 'global_' + ref_name + '.txt'
    elif edist_strategy == 'l':
        ans_file += 'local_' + ref_name + '.txt'
    print "\twriting to '" + ans_file + "'..."
    with open(ans_file, 'w') as outfile:
        outfile.write('>' + ref_name + '\n')
        outfile.write('>INS\n')
        for insert in inserts:
            outfile.write(str(insert[0]) + ',' + str(insert[1]) + '\n')
        outfile.write('>DEL\n')
        for delete in deletes:
            outfile.write(str(delete[0]) + ',' + str(delete[1]) + '\n')
        outfile.write('>SNP\n')
        for snp in snps:
            outfile.write(
                str(snp[0]) + ',' + str(snp[1]) + ',' + str(snp[2]) + '\n')
    return ans_file


def find_variants(ref, ref_name, donor_file, edist_strategy='a'):
    print 'Finding variants...'
    snps = []
    inserts = []
    deletes = []
    donor, donor_name = cm.load_genome(donor_file)
    # loop through 100-length chunks of corresponding ref and donor genome
    chunk_size = 100
    for i in xrange(0, len(ref), chunk_size):
        candidates = []
        for j in xrange(i, min(i + chunk_size, len(ref))):
            if ref[j] != donor[j]:
                candidates.append((ref[j], donor[j], j))
        # NOTE: consider using coefficient of variation/dispersion
        # to determine whether full alignment is needed
        if len(candidates) == 0:    # no variants
            continue
        elif len(candidates) == 1:  # single SNP
            snps.append(candidates[0])
        elif len(candidates) == 2 and candidates[1][2] - candidates[0][2] >= 2:
            # non-consecutive SNPs
            for candidate in candidates:
                snps.append(candidate)
        else:
            # if more than 25% of alleles don't match, don't even bother.
            # We're not likely to find the correct variations
            if len(candidates) > 0.25 * chunk_size:
                continue
            # look for consecutive chunks of mismatches to align
            BUFF_SPACE = 5
            align_start = i
            align_stop = min(i + chunk_size, len(ref))
            # align_start = max(0, candidates[0][2] - BUFF_SPACE)
            # align_stop = min(len(ref), candidates[len(candidates)-1][2] + BUFF_SPACE + 1)

            # try to perform edit-distance alignment on mismatched subsequence
            if edist_strategy == 'a':
                ref_align, donor_align = ed.affine_align(
                    ref[align_start: align_stop],
                    donor[align_start: align_stop]
                )
            elif edist_strategy == 'g':
                ref_align, donor_align = ed.global_align(
                    ref[align_start: align_stop],
                    donor[align_start: align_stop]
                )
            if edist_strategy == 'l':
                ref_align, donor_align = ed.local_align(
                    ref[align_start: align_stop],
                    donor[align_start: align_stop]
                )
            # print 'ref:   ' + ref_align
            # print 'donor: ' + donor_align
            # print
            insert = ''
            delete = ''
            in_insert = False
            in_delete = False
            insert_pos = 0
            delete_pos = 0
            offset = 0
            # ref_align and donor_align equal lengths since we aligned globally
            for j in xrange(len(ref_align)):
                if ref_align[j] == '-':  # insert case
                    if not in_insert:
                        if in_delete:
                            in_delete = False
                            deletes.append((delete, delete_pos))
                            delete = ''
                        in_insert = True
                        insert_pos = i + j - offset
                    insert += donor_align[j]
                    offset += 1
                elif donor_align[j] == '-':  # delete case
                    if not in_delete:
                        if in_insert:
                            in_insert = False
                            inserts.append((insert, insert_pos))
                            insert = ''
                        in_delete = True
                        delete_pos = i + j - offset
                    delete += ref_align[j]
                else:
                    if in_insert:
                        in_insert = False
                        inserts.append((insert, insert_pos))
                        insert = ''
                    elif in_delete:
                        in_delete = False
                        deletes.append((delete, delete_pos))
                        delete = ''
                    if ref_align[j] != donor_align[j]:  # if SNP
                        snps.append(
                            (ref_align[j], donor_align[j], i + j - offset))
        progress = float(i) / float(len(ref))
        cm.print_progress(progress, 'searching')
    print
    return __write_answer_file(ref_name, snps, inserts, deletes, edist_strategy)


def __main():
    args = __parse_args()

    # get the ref and consensus file paths
    # genome_name = cm.get_genome_name(args.genome_name)
    # ref_file = cm.REFS_DIR + cm.REF_PRE + genome_name + '.txt'
    # cons_file = cm.CONS_DIR + cm.CONS_PRE+ genome_name + '.txt'

    ref_file = args.ref_file
    cons_file = args.cons_file

    print 'Loading ref-genome...'
    ref, ref_name = cm.load_genome(ref_file)
    ans_file = find_variants(ref, ref_name, cons_file, args.edit_dist)
    print 'DONE'


if __name__=='__main__':
    __main()