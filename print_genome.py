import argparse
import common as cm


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description="Print a genome by lines of length 100"
    )
    arg_parser.add_argument(
        'genome_file',
        type=str,
        help='The file containing the genome')
    args = arg_parser.parse_args()
    return args


def pretty_print_genome(genome, genome_name):
    print 'Name: ' + genome_name
    line_length = 100
    for i in xrange(0, len(genome), line_length):
        # print
        # print '-' * line_length
        print 'position: ' + str(i)
        print genome[i: i+line_length]
        # print '-' * line_length


def main():
    args = parse_args()
    genome, genome_name = cm.load_genome(args.genome_file)
    pretty_print_genome(genome, genome_name)

if __name__ == '__main__':
    main()