##########################################################
# Name: Ran Moshe Vaknin
# ID: 207761255
# user name: ran.vaknin
# File_name: commands_implement.py
# Description: the file that implement all the different 
# commands that the user chooses with the command lines. 
# all the functions are seperated to different functions.
# this code uses all the functions from 'helper.py'
########################################################## 

##########################################################
#                        imports                         #
##########################################################
import helper
import pickle
import gzip
from classes import Read, Genome
import json


def create_pickle_file(file_name: str, data: dict) -> None:
    '''creates pickle compressed file with a given data and given name and returns none
    in case there is no file name'''
    try:
        with gzip.open(file_name, "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    except TypeError:
        print(f'An output file name is missing!')
    return None


def load_pickle_file(file_name: str) -> dict[str, dict[str, dict] | list[Genome]] | None:
    '''returns the data from a given pickle file and returns none in case of it's not exist'''
    try:
        with gzip.open(file_name, 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError:
        print(f'{file_name} does not exist.')
    return None


def reference_command_implementation(genomes_file: str, KMER_SIZE: str, output_kdb_file: str,
                                     ext_filter_similar: bool, ext_similarity_threshold: str) -> None:
    '''Implement the command reference. Creates a pickle file with the data from the function 'genomes_summerize'
    from helper.'''
    try:
        if output_kdb_file==None:
            raise FileNotFoundError
        KMER_SIZE = int(KMER_SIZE)
    except FileNotFoundError:
        print (f'An output file name is missing!')
    except TypeError:
        print('kmer_size is missing! Use the flag of -k')
    else:
        ext_similarity_threshold = float(ext_similarity_threshold)
        all_genomes_data = helper.genomes_summerize(
            genomes_file, KMER_SIZE, ext_filter_similar, ext_similarity_threshold)
        if all_genomes_data == None:
            # there's an error with the genomes file
            return None
        create_pickle_file(output_kdb_file, all_genomes_data)
    return None


def dumpref_command_implementation(genomes_file: str | None, KMER_SIZE: str | None, reference_file: str | None,
                                   ext_filter_similar: bool, ext_similarity_threshold: str) -> None:
    '''Implement the command dumpref. Seperate between two cases:
    1. getting only referrence file 
    2. FASTA file and kmer size. 
    Implement the highly similar genomes filterring extention if needed.'''
    all_genomes_data = None
    ext_similarity_threshold = float(ext_similarity_threshold)
    try:
        if reference_file != None:
            # case1
            all_genomes_data = load_pickle_file(reference_file)
        elif KMER_SIZE != None and genomes_file != None:
            # case2
            KMER_SIZE = int(KMER_SIZE)
            all_genomes_data = helper.genomes_summerize(
                genomes_file, KMER_SIZE, ext_filter_similar, ext_similarity_threshold)
        else:
            # if the input is not as needed
            raise ValueError
    except ValueError:
        print(f'Invalid input!')
    else:
        if all_genomes_data != None:
            similarity_results = all_genomes_data['similarity_results']
            references_data = all_genomes_data['references_data']
            kmers_details = references_data['kmers_data']
            genomes_list = references_data['genomes_data']
            genomes_summary = dict()
            for genome in genomes_list:
                genomes_summary[genome.get_name(
                )] = genome.get_genome_summary()
            if ext_filter_similar:
                # prints the extention data if needed
                references_summary = {'Kmers': kmers_details, 'Summary': genomes_summary,'Similarity': similarity_results}
                print(json.dumps(references_summary))
            else:
                references_summary = {
                    'Kmers': kmers_details, 'Summary': genomes_summary}
                print(json.dumps(references_summary))                
    return None


def reference_data_creator(reference_file: str | None, genomes_file: str | None,
                           KMER_SIZE: str | None) -> tuple[list[Genome] | None, dict[str, dict] | None, int | None]:
    '''Returns all the data from a reference file or the 'genomes_summerize' function depends on the 
    given data'''
    genomes_list = None
    kmers_details = None
    if reference_file != None:
        all_genomes_data = load_pickle_file(reference_file)
        if all_genomes_data != None:
            references_data = all_genomes_data['references_data']
            kmers_details = references_data['kmers_data']
            genomes_list = references_data['genomes_data']
            KMER_SIZE = len(list(kmers_details.keys())[0])
    elif genomes_file != None and KMER_SIZE != None:
        KMER_SIZE = int(KMER_SIZE)
        all_genomes_data = helper.genomes_summerize(genomes_file, KMER_SIZE)
        if all_genomes_data != None:
            references_data = all_genomes_data['references_data']
            kmers_details = references_data['kmers_data']
            genomes_list = references_data['genomes_data']
    else:
        raise ValueError
    return genomes_list, kmers_details, KMER_SIZE


def reads_data_creator(reads_list: list[Read], kmers_details: dict[str, dict], genomes_list: list[Genome],
                       UNIQUE_THRESHOLD: str, AMBIGUOUS_THRESHOLD: str, min_kmer_quality: str | None,
                       max_genomes: str | None, ext_quality_summary: dict[str, int]) -> dict[str, dict[str, int] | list[Genome], dict[str, int]]:
    '''Returns all the reads data after derminate it by step 1 and step 2 (if ambiguous threshold given) 
    functions from helper. Returns the reads statistics, list of all genomes and ext_quality_summary if 
    one of the quality extention argument is given.'''
    if min_kmer_quality != None:
        min_kmer_quality = int(min_kmer_quality)
    if max_genomes != None:
        max_genomes = int(max_genomes)
    filtered_sorted_reads = helper.reads_kmers_specifier(
        kmers_details, reads_list, min_kmer_quality, max_genomes)
    # A dictionary with all the needed data from 'reads_kme_specifier' function from helper
    ext_quality_summary['filtered_quality_kmers'] = filtered_sorted_reads['filtered_quality_kmers']
    ext_quality_summary['filtered_hr_genomes'] = filtered_sorted_reads['filtered_hr_genomes']
    reads_list = filtered_sorted_reads['reads_list']
    reads_list, reads_stats, genomes_list = helper.step1_determination_definer(
        reads_list, UNIQUE_THRESHOLD, genomes_list)
    if AMBIGUOUS_THRESHOLD.isdigit():
        AMBIGUOUS_THRESHOLD = int(AMBIGUOUS_THRESHOLD)
        _, reads_stats, genomes_list = helper.step2_determination_definer(
            reads_list, reads_stats, AMBIGUOUS_THRESHOLD, genomes_list)
    reads_data = {'reads_stats': reads_stats, 'genomes_list': genomes_list,
                  'ext_quality_summary': ext_quality_summary}
    return reads_data


def align_command_implementation(reference_file: str | None, reads_file: str | None, UNIQUE_THRESHOLD: str,
                                 AMBIGUOUS_THRESHOLD: str, genomes_file: str | None, KMER_SIZE: str | None, align_file: str | None,
                                 min_read_quality: str | None, min_kmer_quality: str | None, max_genomes: str | None) -> None:
    '''Implements the 'align' command and creates a pickle compressed file with all the needed data of the 
    reads and all the genomes. Implement the quality extentiokn if one of the arguments is given'''
    try:
        if align_file==None:
            raise FileNotFoundError
        genomes_list, kmers_details, KMER_SIZE = reference_data_creator(
            reference_file, genomes_file, KMER_SIZE)
    except FileNotFoundError:
        print(f'An output file name is missing!')    
    except ValueError:
        # in case the given input is not as needed
        print(f'Invalid input!')
    else:
        if genomes_list == None:
            return None
        if min_read_quality != None:
            min_read_quality = int(min_read_quality)
        reads_data_n_ext_quality = helper.fastq_reader(
            reads_file, KMER_SIZE, min_read_quality)
        if reads_data_n_ext_quality == None:
            # in case there is no FASTQ valid file
            return None
        reads_list = reads_data_n_ext_quality['reads_list']
        ext_quality_summary = reads_data_n_ext_quality['ext_quality_summary']
        UNIQUE_THRESHOLD = int(UNIQUE_THRESHOLD)
        reads_data = reads_data_creator(
            reads_list, kmers_details, genomes_list, UNIQUE_THRESHOLD, AMBIGUOUS_THRESHOLD,
            min_kmer_quality, max_genomes, ext_quality_summary)
        create_pickle_file(align_file, reads_data)
    return None


def dumpalign_command_implementation(align_file: str | None, reference_file: str | None, reads_file: str | None, UNIQUE_THRESHOLD: str | None,
                                     AMBIGUOUS_THRESHOLD: str | None, genome_file: str | None, KMER_SIZE: str | None,
                                     min_read_quality: str | None, min_kmer_quality: str | None, max_genomes: str | None) -> None:
    '''Implement the 'dumpalign' command by the given data and implement the quality extentiokn if one of the arguments is given'''
    try:
        if align_file != None:
            reads_data = load_pickle_file(align_file)
        elif reads_file != None:
            genomes_list, kmers_details, KMER_SIZE = reference_data_creator(
                reference_file, genome_file, KMER_SIZE)
            if genomes_list == None:
                return None
            if min_read_quality != None:
                min_read_quality = int(min_read_quality)
            reads_data_n_ext_quality = helper.fastq_reader(
                reads_file, KMER_SIZE, min_read_quality)
            UNIQUE_THRESHOLD = int(UNIQUE_THRESHOLD)
            if reads_data_n_ext_quality == None:
                return None
            reads_list = reads_data_n_ext_quality['reads_list']
            ext_quality_summary = reads_data_n_ext_quality['ext_quality_summary']
            reads_data = reads_data_creator(reads_list, kmers_details, genomes_list,
                                            UNIQUE_THRESHOLD, AMBIGUOUS_THRESHOLD, min_kmer_quality, max_genomes, ext_quality_summary)
        else:
            raise ValueError
    except ValueError:
        print(f'Invalid input!')
    else:
        if reads_data == None:
            return None
        ext_quality_summary = reads_data['ext_quality_summary']
        if min_read_quality == None:
            ext_quality_summary.pop('filtered_quality_reads')
        if min_kmer_quality == None:
            ext_quality_summary.pop('filtered_quality_kmers')
        if max_genomes == None:
            ext_quality_summary.pop('filtered_hr_genomes')
        genomes_list = reads_data['genomes_list']
        reads_stats = reads_data['reads_stats']
        genome_mapping_summary = {}
        for genome in genomes_list:
            genome_mapping_summary[genome.get_name(
            )] = genome.get_mapping_summary()
        final_reads_data = {'Statistics': {**reads_stats, **ext_quality_summary},
                            'Summary': genome_mapping_summary}
        print(json.dumps(final_reads_data))
