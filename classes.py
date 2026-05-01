##########################################################
# Name: Ran Moshe Vaknin
# ID: 207761255
# user name: ran.vaknin
# File_name: classes.py
# Description: This code contains all the basic objects 
# for the whole program. 
# Classes names- KmerReference, Genome,Read
########################################################## 

##########################################################
#                        imports                         #
##########################################################
import numpy as np


class KmerReference:
    '''A class that represent for each kmer in the 
    reference file the data, such as the kmer name and the 
    indices that the kmer appears in each genome .'''

    def __init__(self, name: str):
        self.__name = name
        self.__genomes_appearances: dict = dict()

    def __repr__(self) -> str:
        return str(self.__genomes_appearances)

    ###############
    #   setters   #
    ###############

    def add_reference(self, gene_name: str, appearances: list[int]) -> None:
        '''Add to the kmer all the indices that it appears in genome'''
        self.__genomes_appearances[gene_name] = appearances

    ###############
    #   getters   #
    ###############

    def get_name(self) -> str:
        return self.__name

    def get_appearances(self) -> dict[str, list[int]]:
        return self.__genomes_appearances


class Genome:
    '''A class that represent each genome from the reference file and 
    the relevant data for the genome'''

    def __init__(self, name: str, total_length: int):
        self.__name = name
        self.__kmers_list: list[str] = []
        self.__genome_summary = {
            'total_bases': total_length, 'unique_kmers': 0, 'multi_mapping_kmers': 0}
        self.__mapping_summary = {'unique_reads': 0 , 'ambiguous_reads': 0}

    def __eq__(self, value: str) -> bool:
        return self.__name == value

    ###############
    #   setters   #
    ###############
    def create_kmers(self, all_genome: str, k: int) -> None:
        '''Creates the kmers list of the genome'''
        genome_array = np.array(list(all_genome))
        self.__kmers_list = [''.join(genome_array[i:i+k])
                             for i in range(len(genome_array)-k+1)]

    def add_unique_read(self) -> None:
        self.__mapping_summary['unique_reads'] += 1

    def add_ambiguous_read(self) -> None:
        self.__mapping_summary['ambiguous_reads'] += 1

    def reduce_unique_read(self) -> None:
        self.__mapping_summary['unique_reads'] -= 1

    def add_unique_kmer(self) -> None:
        self.__genome_summary['unique_kmers'] += 1

    def add_multi_mapping_kmer(self) -> None:
        self.__genome_summary['multi_mapping_kmers'] += 1

    ###############
    #   getters   #
    ###############
    def get_kmers_list(self) -> list[str]:
        return self.__kmers_list

    def get_name(self) -> str:
        return self.__name

    def get_genome_summary(self) -> dict[str, int]:
        return self.__genome_summary

    def get_mapping_summary(self) -> dict[str, int]:
        return self.__mapping_summary


class Read:
    '''A calss that represents each read from the fastq file and the 
    relevant data for this read'''

    def __init__(self, name: str):
        self.__name = name
        self.__kmers_list: list[str] = []
        self.__quality = None
        self.__specific_kmers = dict()
        self.__unspecific_kmers = dict()
        self.__determination = ''

    ###############
    #   setters   #
    ###############
    def create_kmers(self, all_read: str, k: int) -> None:
        '''creates a list of the read's kmers'''
        read_array = np.array(list(all_read))
        self.__kmers_list = [''.join(read_array[i:i+k])
                             for i in range(len(read_array)-k+1)]

    def add_specific_kmer(self, kmer: str, genomes_appearances: dict[str, list[int]]) -> None:
        self.__specific_kmers[kmer] = genomes_appearances

    def add_unspecific_kmer(self, kmer: str, genomes_appearances: dict[str, list[int]]) -> None:
        self.__unspecific_kmers[kmer] = genomes_appearances

    def define_quality(self, quality: str) -> None:
        self.__quality = quality

    def define_determination(self, determination: str) -> None:
        self.__determination = determination

    def remove_kmer(self, kmer: str) -> None:
        '''removes a kmer from the kmers list caused by filtering extention'''
        self.__kmers_list.remove(kmer)

    ###############
    #   getters   #
    ###############
    def create_quality_by_kmers(self) -> list[str]:
        if len(self.__kmers_list) == 0:
            kmer_size = 0
        else:
            kmer_size = len(self.__kmers_list[0])
        quality_array = np.array(list(self.__quality))
        quality_by_kmers = [''.join(quality_array[i:i+kmer_size])
                            for i in range(len(quality_array)-kmer_size+1)]
        return quality_by_kmers

    def get_specific_kmers(self) -> dict[str, dict[str, int]]:
        return self.__specific_kmers

    def get_unspecific_kmers(self) -> dict[str, dict[str, int]]:
        return self.__unspecific_kmers

    def get_name(self) -> str:
        return self.__name

    def get_determination(self) -> str:
        return self.__determination

    def get_kmers_list(self) -> list[str]:
        return self.__kmers_list

    def get_quality(self) -> str:
        return self.__quality
