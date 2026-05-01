##########################################################
# Name: Ran Moshe Vaknin
# ID: 207761255
# user name: ran.vaknin
# File_name: helper.py
# Description: The code that do all the different needed
# functions in this program. all the data created through 
# this code by using its different functions.
########################################################## 

##########################################################
#                        imports                         #
##########################################################
from classes import Genome, KmerReference, Read
import re
from collections import defaultdict


def fasta_reader(file_path: str) -> dict[str, str]:
    '''reads a file in fasta format and returns a dictionary 
    in the format {gene_name: whole_sequence}'''
    sequence = dict()
    current_header = None
    fasta_file = open(file_path)
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            current_header = line[1:]
            sequence[current_header] = ''
        else:
            sequence[current_header] += line
    return sequence


def fastq_reader(file_path: str, KMER_SIZE: int, min_read_quality: int | None) -> dict[str, list[Read] | dict[str, int]] | None:
    '''Reads a FASTQ file and creates a list of Read objects.
    for each read it will define it's name and the quality of it and it will 
    create a list of all the kmers of the read as kmers_list in the object.'''
    ext_quality_summary = {'filtered_quality_reads': 0}
    # for the extention EXTQUALITY
    try:
        fastq_file = iter(open(file_path))
    except FileNotFoundError or TypeError:
        print(f'Enter a valid FASTQ file.')
        return None
    # stops the whole progran in case of the file name doesn't exist
    else:
        reads_list = []
        current_read = None
        default = 'The end of the file'
        line = next(fastq_file, default)
        while line is not default:
            line = line.strip()
            current_read = Read(line[1:])
            reads_list.append(current_read)
            line = next(fastq_file, default)
            # the read's sequence
            line = line.strip()
            current_read.create_kmers(line, KMER_SIZE)
            line = next(fastq_file, default)
            # the '+' line
            line = next(fastq_file, default)
            # the read's quality line
            line = line.strip()
            current_read.define_quality(line)
            if min_read_quality != None:
                # in case of minimum read's quality filtering extension is active it
                # will check the quality of the whole read and will ignore from it if
                # it lower than the threshold
                read_quality = sum(
                    ord(i) - 33 for i in current_read.get_quality())/len(current_read.get_quality())
                if read_quality < min_read_quality:
                    reads_list.remove(current_read)
                    ext_quality_summary['filtered_quality_reads'] += 1
            line = next(fastq_file, default)
        return {'reads_list': reads_list, 'ext_quality_summary': ext_quality_summary}


def genomes_list_creator(genomes: dict[str, str], k: int) -> list[Genome]:
    '''Creates a list of Genome objects from the FASTA file data'''
    all_genes = []
    for genome in genomes:
        current_gene = Genome(genome, len(genomes[genome]))
        all_genes.append(current_gene)
        current_gene.create_kmers(genomes[genome], k)
    return all_genes


def kmers_references(file_path: str, KMER_SIZE: int) -> tuple[dict[str, KmerReference], list[Genome]]:
    """Creates a dictionary mapping each k-mer to its reference information.
    Returns a dictionary of al the kmers references and list of all the genomes"""
    all_sequences = fasta_reader(file_path)
    all_genomes = genomes_list_creator(all_sequences, KMER_SIZE)
    kmer_positions: defaultdict[defaultdict[str, list]
                                ] = defaultdict(lambda: defaultdict(list))
    # creates this type of dict {kmer: {gene_name: [indices]}}
    for genome in all_genomes:
        kmers_list = genome.get_kmers_list()
        gene_name = genome.get_name()
        for idx, kmer in enumerate(kmers_list):
            if bool(re.match('^[AGCT]+$', kmer)):
                kmer_positions[kmer][gene_name].append(idx)
                # adds a value of the type {kmer_name: {genome_name: [idx]}}
    all_kmer_references = {
        kmer: KmerReference(kmer) for kmer in kmer_positions
    }
    # creates a dictionary like {kmer_name: KmerRefrence}
    for kmer, genome_refs in kmer_positions.items():
        for genome, positions in genome_refs.items():
            all_kmer_references[kmer].add_reference(genome, positions)
            # adds a dict {genome_name: [positions]} to each kmer
    return all_kmer_references, all_genomes


def update_kmer_details(kmers_details: dict[str, dict[str, list[int]]], genome_name: str) -> dict[str, dict]:
    '''Removes a genomoe frome the kmer details dictionary that looks like: {kmer:{genome_name:[indices]}}'''
    for kmer_appearances in kmers_details:
        if genome_name in kmers_details[kmer_appearances]:
            kmers_details[kmer_appearances].pop(genome_name)
    return kmers_details


def similarity_filtering(genomes_list: list[Genome], kmers_details: dict[str, dict[str, list[int]]],
                         ext_similarity_threshold: float) -> tuple[list[Genome], dict[str, dict[str, list[int]]], dict[str, str | int]]:
    '''Filters the kmer_details dictionary by removing genomes by the algorithm of 
    filter_similar extention
    return: list of Genome objects, kmers dictionatry like {kmer:{genome_name:[indices]}}, similarity results'''
    similarity_results = dict()
    sorted_genomes_list = sorted(genomes_list, key=lambda genome: (genome.get_genome_summary()['unique_kmers'],
                                                                   len(genome.get_kmers_list()), genome.get_genome_summary()['total_bases'], genome.get_name()))
    # creates a copy of the genomes list and sorts it by this levels:
    # specific kmers count -> total kmers count -> genome sequence length -> genome name
    for genome_1_idx, genome_1 in enumerate(sorted_genomes_list):
        genome_1_kmers = set(genome_1.get_kmers_list())
        genome_1_removed = False
        for genome_2_idx, genome_2 in enumerate(sorted_genomes_list):
            if genome_1_idx == genome_2_idx:
                continue
            # in case genome_2 is genome_1 itself
            genome_2_kmers = set(genome_2.get_kmers_list())
            similarity_score = len(genome_1_kmers.intersection(
                genome_2_kmers))/min(len(genome_1_kmers), len(genome_2_kmers))
            # calculates the similarity score f the current genome with each other genome
            # by the common kmers of the two genomes count divided by the lower size of
            # the total kmers count
            if similarity_score > ext_similarity_threshold:
                genomes_list.remove(genome_1)
                similarity_results[genome_1.get_name()] = {'kept': 'no', 'unique_kmers': genome_1.get_genome_summary()['unique_kmers'],
                                                           'total_kmers': len(genome_1.get_kmers_list()), 'genome_length': genome_1.get_genome_summary()['total_bases'],
                                                           'similar_to': genome_2.get_name(), 'similarity_score': similarity_score}
                kmers_details = update_kmer_details(
                    kmers_details, genome_1.get_name())
                # removes the genome in case of the similarity score is lower than the similarity threshold
                genome_1_removed = True
                break
                # continue to the next genome for examination
        if genome_1_removed:
            continue
        similarity_results[genome_1.get_name()] = {'kept': 'yes', 'unique_kmers': genome_1.get_genome_summary()['unique_kmers'],
                                                   'total_kmers': len(genome_1.get_kmers_list()), 'genome_length': genome_1.get_genome_summary()['total_bases'],
                                                   'similar_to': 'NA', 'similarity_score': 'NA'}
    return genomes_list, kmers_details, similarity_results


def genomes_summerize(file_path: str, KMER_SIZE: int, ext_filter_similar: bool = None, ext_similarity_threshold: int = None) -> dict[str, dict | None] | None:
    '''Summerize all the needed information of the reference file included the information for the similarity filter extention.
     return a dictionary with the refrence database (such as the appearances of the kmers in genomes and the genomes list)
      and the relevant information for the similarity filter extention.'''
    similarity_results = None
    try:
        all_kmers_references, genomes_list = kmers_references(
            file_path, KMER_SIZE)
    except FileNotFoundError or TypeError:
        print('Enter a valid file in FASTA format.')
        # in case the FASTA file doesn't exist the program will stop.
    else:
        genomes_summery = dict()
        kmers_details = dict()
        for genome in genomes_list:
            genomes_summery[genome.get_name()] = genome.get_genome_summary()
            # creates dictionary like {genome_name: genome_summary}
        for kmer in all_kmers_references.values():
            kmers_details[kmer.get_name()] = kmer.get_appearances()
            # creates a dictionary like {kmer:{genome_name:[indices]}}
            if len(kmer.get_appearances()) == 1:
                for genome in kmer.get_appearances():
                    genome_idx = genomes_list.index(genome)
                    genomes_list[genome_idx].add_unique_kmer()
                    # in case the kmer appears only in one genome
            else:
                for genome in kmer.get_appearances():
                    genome_idx = genomes_list.index(genome)
                    genomes_list[genome_idx].add_multi_mapping_kmer()
        if ext_filter_similar:
            genomes_list, kmers_details, similarity_results = similarity_filtering(
                genomes_list, kmers_details, ext_similarity_threshold)
        return {'references_data': {'kmers_data': kmers_details, 'genomes_data': genomes_list}, 'similarity_results': similarity_results}


def reads_kmers_specifier(kmers_refrences: dict[str, dict], reads: list[Read], min_kmer_quality: int | None, max_genomes: int | None) -> dict[str, list[Read] | int]:
    '''Adds to each read a specific kmers and unspecific kmers by the genomes they appear in
    and fiter the kmers in each read by kmers minimum quality and max mapped genomes extentions.
    returns a dictionary with list of all the reads as Read objects and the amount of kmers that 
    filtered by each extention (if it's not active it will be 0)'''
    reads_data = {'reads_list': reads,
                  'filtered_quality_kmers': 0, 'filtered_hr_genomes': 0}
    for read in reads:
        kmers_list_copy = read.get_kmers_list()
        if min_kmer_quality != None:
            quality_by_kmers = read.create_quality_by_kmers()
        for kmer_idx, kmer in enumerate(kmers_list_copy):
            if min_kmer_quality != None:
                kmer_quality = sum(
                    ord(i) - 33 for i in quality_by_kmers[kmer_idx])/len(kmer)
                if kmer_quality < min_kmer_quality:
                    read.remove_kmer(kmer)
                    reads_data['filtered_quality_kmers'] += 1
                    continue
                # filters the kmers of the current read by the kmers minimum quality extention
            if bool(re.match('^[AGCT]+$', kmer)):
                # checks if there are no letters except 'AGCT' in the kmer. if so it ignores from the kmer.
                if kmer in kmers_refrences:
                    kmer_details = dict()
                    for gene in kmers_refrences[kmer]:
                        kmer_details[gene] = len(
                            kmers_refrences[kmer][gene])
                        # adds the amount that the kmer appears in the genome
                    if len(kmer_details) == 1:
                        if max_genomes != None and max_genomes <= 1:
                            read.remove_kmer(kmer)
                            reads_data['filtered_hr_genomes'] += 1
                            continue
                            # remove kmer by the max mapped genomes extention if it's activate
                        read.add_specific_kmer(kmer, kmer_details)
                    else:
                        if max_genomes != None and len(kmer_details) > max_genomes:
                            read.remove_kmer(kmer)
                            reads_data['filtered_hr_genomes'] += 1
                            continue
                            # remove kmer by the max mapped genomes extention if it's activate
                        read.add_unspecific_kmer(kmer, kmer_details)
    return reads_data


def specific_kmers_one_genome(specific_kmers: dict[str, dict]) -> tuple[bool, str | None]:
    '''checks if all the specific kmers in a read are mapping only one genome.
    if so it returns the genome name '''
    genomes = {list(genomes.keys())[0] for genomes in specific_kmers.values()}
    return len(genomes) == 1, genomes.pop() if len(genomes) == 1 else None


def specific_genomes_sorted_list(specific_kmers: dict[str, dict]) -> list[tuple[str, int]]:
    """Creates sorted list of genomes based on k-mer counts, gets a dictionary of 
    specific kmers ans their genomes mapping and returns list of tuples like
    (genome_name, k-mer_count) sorted by count."""
    genome_counts = defaultdict(lambda: 0)
    for genomes in specific_kmers.values():
        for genome, count in genomes.items():
            genome_counts[genome] += count
    # Count specific k-mers per genome
    sorted_genomes = sorted(genome_counts.items(),
                            key=lambda x: x[1], reverse=True)
    # Get the top two genomes
    return sorted_genomes


def step1_determination_definer(reads: list[Read], UNIQUE_THRESHOLD: int,
                                genomes_list: list[Genome]) -> tuple[list[Read], dict[str, int], list[Genome]]:
    '''Determinates each read by its specific kmers count. Determines by the given unique threshold.
    Returns list of all reads as Read objects, statistics of all reads determinations and all genomes 
    list as Genome objects.'''
    reads_stats = {'unique_mapped_reads': 0,
                   'ambiguous_mapped_reads': 0, 'unmapped_reads': 0}
    for read in reads:
        if len(read.get_specific_kmers()) == 0:
            read.define_determination('UNMAPPED')
            reads_stats['unmapped_reads'] += 1
            continue
        specific_kmers = read.get_specific_kmers()
        is_single_genome, genome_name = specific_kmers_one_genome(
            specific_kmers)
        # checks if the read's kmers map only one genome
        if is_single_genome == True:
            read.define_quality(f'UNIQUELY MAPPED to {genome_name}')
            reads_stats['unique_mapped_reads'] += 1
            genome_idx = genomes_list.index(genome_name)
            genomes_list[genome_idx].add_unique_read()
            continue
        sorted_genomes_list = specific_genomes_sorted_list(specific_kmers)
        if sorted_genomes_list[0][1] - sorted_genomes_list[1][1] >= UNIQUE_THRESHOLD:
            # if the difference between the count of specific kmers of the highest
            # two genomes is higher or equal to the given unique threshold, read will
            # be clasify as uniquely mapped and otherwise as ambiguously mapped
            read.define_determination(
                f'UNIQUELY MAPPED to {sorted_genomes_list[0][0]}')
            genome_idx = genomes_list.index(sorted_genomes_list[0][0])
            genomes_list[genome_idx].add_unique_read()
            reads_stats['unique_mapped_reads'] += 1
            continue
        else:
            specific_mapping_genomes = []
            for genome in sorted_genomes_list:
                specific_mapping_genomes.append(genome[0])
            read.define_determination(
                f'AMBIGUOUSLY MAPPED')
            for genome in specific_mapping_genomes:
                genome_idx = genomes_list.index(genome)
                genomes_list[genome_idx].add_ambiguous_read()
            reads_stats['ambiguous_mapped_reads'] += 1
    return reads, reads_stats, genomes_list


def only_uniquely_reads(reads: list[Read]) -> list[Read]:
    '''Return a list of all the uniquely mapped reads'''
    filtered_list = []
    for read in reads:
        if 'UNIQUELY MAPPED' in read.get_determination():
            filtered_list.append(read)
    return filtered_list


def read_genomes_databse(read: Read) -> dict[str, dict[str, int]]:
    '''Returns a dictionary like {genome_name:
    {specific:(count of specificaly map in the read), total:(count of total maps in the read)}
    }'''
    genomes_db = defaultdict(lambda: {'specific': 0, 'total': 0})
    for genomes in read.get_specific_kmers().values():
        for genome in genomes.keys():
            genomes_db[genome]['specific'] += 1
            genomes_db[genome]['total'] += 1
    for genomes in read.get_unspecific_kmers().values():
        for genome, _ in genomes.items():
            genomes_db[genome]['total'] += 1
    return genomes_db


def step2_determination_definer(reads: list[Read], reads_stats: dict[str, int], AMBIGUOUS_THRESHOLD: int,
                                genomes_list: list[Genome]) -> tuple[list[Read], dict[str, int], list[Genome]]:
    '''Determinates each uniquely mapped read by the algoritm of the ambiguous threshold.
    Returns list of all reads as Read objects, statistics of all reads determinations and all genomes 
    list as Genome objects.'''
    reads_filtered = only_uniquely_reads(reads)
    for read in reads_filtered:
        genomes_mapping = read_genomes_databse(read)
        max_count_genome_values = max(
            genomes_mapping.values(), key=lambda value: value['total'])
        # gets the values of the genome with the highest count of total kmers
        map_count_genome_values = max(
            genomes_mapping.values(), key=lambda value: value['specific'])
        # gets the values of the genome with the highest count of specific kmers
        if max_count_genome_values['total']-map_count_genome_values['total'] > AMBIGUOUS_THRESHOLD:
            # in case the difference between the read with the highest kmer counts and the
            # kmer counts of the read with the highest uniquely mapped kmers count is
            # higher from the given ambiguous threshold
            map_count_value = map_count_genome_values['total']
            ambiguous_genomes_list = [genome for genome, values in genomes_mapping.items(
            ) if values['total'] >= map_count_value]
            uniquely_mapped_genome = read.get_determination().replace('UNIQUELY MAPPED to ', '')
            genome_idx = ambiguous_genomes_list.index(uniquely_mapped_genome)
            genomes_list[genome_idx].reduce_unique_read()
            read.define_determination(f'AMBIGUOUSLY MAPPED')
            for ambiguous_genome in ambiguous_genomes_list:
                genome_idx = genomes_list.index(ambiguous_genome)
                genomes_list[genome_idx].add_ambiguous_read()
            reads_stats['unique_mapped_reads'] -= 1
            reads_stats['ambiguous_mapped_reads'] += 1
    return reads, reads_stats, genomes_list
