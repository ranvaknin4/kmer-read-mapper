#!/usr/bin/env python3
##########################################################
# Name: Ran Moshe Vaknin
# ID: 207761255
# user name: ran.vaknin
# File_name: main.py
# Description: The main part that define all the possible
# functions and use the functions in 'commands_implement.py'
# to run the program with the relevant argument for each 
# command.
########################################################## 

##########################################################
#                        imports                         #
##########################################################
import argparse
import commands_implement as cmd_imp


def readargs(args=None):
    parser = argparse.ArgumentParser(
        prog='Biosequence project',
    )
    # General arguments
    parser.add_argument('-t', '--task',
                        help="task", choices=['reference', 'dumpref', 'align', 'dumpalign'],
                        required=True
                        )
    parser.add_argument('-g', '--genomefile',
                        help="Genome fasta file (multiple records)",
                        )
    parser.add_argument('-r', '--referencefile',
                        help="kdb file. Can be either input or name for output file",
                        )
    parser.add_argument('-k', '--kmer-size', type=is_positive_number,
                        help="length of kmers",
                        )
    # align+dumpalign
    parser.add_argument('-a', '--alignfile',
                        help="aln file. Can be either input or name for output file",
                        )
    parser.add_argument('--reads',
                        help="fastq reads file",
                        )
    parser.add_argument('-m', '--unique-threshold', type=is_positive_number, default='1',
                        help="unique k-mer threshold",
                        )
    parser.add_argument('-p', '--ambiguous-threhold', default='1',
                        help="ambiguous k-mer threshold",
                        )
    # align+quality
    parser.add_argument('--min-read-quality',
                        type=int,
                        )
    parser.add_argument('--min-kmer-quality',
                        type=int,
                        )
    parser.add_argument('--max-genomes',
                        type=int,
                        )
    # similarity
    parser.add_argument('--filter-similar',
                        action='store_true'
                        )
    parser.add_argument('--similarity-threshold',
                        type=is_interval,
                        default=0.95,
                        )

    return parser.parse_args(args)


def is_positive_number(value: str) -> str:
    '''Cheks if argumants that should be integers are and raise an error 
    if they are not.'''
    try:
        int_value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f'{value} is invalid. Must be a non-negative integer.')
    if int_value < 0:
        raise argparse.ArgumentTypeError(
            f'{value} is invalid. Must be a non-negative integer.')
    return value


def is_interval(value: str) -> str:
    '''checks if a given value is a number between 0 to 1 and raise error if not '''
    try:
        float_value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f'{value} is invalid. Must be an interval(0<=value<=1).')
    if float_value < 0 or float_value > 1:
        raise argparse.ArgumentTypeError(
            f'{value} is invalid. Must be an interval(0<=value<=1).')
    return value


def main() -> None:
    args = readargs()
    if args.task == 'reference':
        cmd_imp.reference_command_implementation(
            args.genomefile, args.kmer_size, args.referencefile, args.filter_similar, args.similarity_threshold)
    elif args.task == 'dumpref':
        cmd_imp.dumpref_command_implementation(
            args.genomefile, args.kmer_size, args.referencefile, args.filter_similar, args.similarity_threshold)
    elif args.task == 'align':
        cmd_imp.align_command_implementation(args.referencefile, args.reads, args.unique_threshold,
                                       args.ambiguous_threhold, args.genomefile, args.kmer_size, args.alignfile,
                                       args.min_read_quality, args.min_kmer_quality, args.max_genomes)
    elif args.task == 'dumpalign':
        cmd_imp.dumpalign_command_implementation(args.alignfile, args.referencefile, args.reads,
                                           args.unique_threshold, args.ambiguous_threhold, args.genomefile,
                                           args.kmer_size, args.min_read_quality, args.min_kmer_quality, args.max_genomes)


if __name__ == "__main__":
    main()
