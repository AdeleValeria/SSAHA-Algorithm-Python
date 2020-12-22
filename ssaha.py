"""
Adele Valeria
Write a SSAHA Algorithm Using Python
Biomedical Informatics 3 Mini-Project

The algorithm was introduced by Zemin Ning et al.
SSAHA: A Fast Search Method for Large DNA Databases
www.genome.org/cgi/doi/10.1101/gr.194201

Programmed using Python 3.8.5

Some inspirations come from:
github.com/Lygeri-Sakellaridi/Custom-Implementations-of-Bioinformatics-Algorithms/blob/e2ba8063bb916cae38c033941e38fa2e381a0d02/ssaha.py

"""

# Import Packages
"""
tabulate 0.8.7 : print output table
requests, sys : interact with Ensembl REST API
optparse 1.5.3: command line options

"""
from tabulate import tabulate
import requests, sys
import optparse

# to track memory usage and processing time
#import os, psutil
#import time
#process = psutil.Process(os.getpid())
#tic = time.clock()

"""
#------RANDOM QUERY SEQUENCES GENERATOR------
def random_dna_sequence(length):
    import numpy as np
    bases = ('A', 'C', 'G', 'T')
    return ''.join(np.random.choice(bases, length))

with open("query.txt","w") as query_text:
    for i in range(20):
        query_text.write(">query #" + str(i+1))
        query_text.write("\n")
        query_text.write(random_dna_sequence(500))
        query_text.write("\n")
        query_text.write("\n")
"""

def reverseComplement(query):
    """
    query : query sequence in "reverse" direction
    
    """
    complement = {'A': 'T', 
                  'C': 'G', 
                  'G': 'C', 
                  'T': 'A'}
    return "".join(complement[i] for i in reversed(query))

def parseFasta(file):
    """
    file : user's input file to be parsed (ref genome and or query)

    Returns
    -------
    label : a list of headers in the file (e.g. [1 dna:chromosome, 2 dna:chromosome, ...])
    all_seq : a list of sequences that have been combined (e.g. [seq1, seq2, ...])
    
    When a new header is encountered for the 2nd time and so on, 
    join sequences in temp list into one string and append it
    to all_seq. Then, empty temp, so it can temporarily retain next sequences.
    
    The file does not end with a header, so the last sequences in temp are
    joined after the loop is completed.
    """
    with open(file, mode = "r") as f:
        label = []
        all_seq = []
        temp = []
        # iterate the file line by line
        for line in f:
            # remove trailing whitespace 
            line = line.rstrip()
            # line beginning with ">" is a header
            if line.startswith(">"):
                # splice[start:stop:step]
                label.append(line[1:])
                if temp != []:
                    # sequences in temp list are separated by comma
                    # so we join them into one string
                    all_seq.append("".join(temp))
                    temp = []
            # any other lines are sequences, so we append each line into a temporary list
            else:
                temp.append(line) 
        chr_seq = "".join(temp)
        all_seq.append(chr_seq)
        temp = []
    return label, all_seq


def hashTable(all_seq, k):
    """
    all_seq : a list of parsed sequences from reference genome
    k : length of kmer, default value is 5 if not specified by user

    Returns
    -------
    
    hash_dict : a dictionary containing non-overlapping kmers (keys) and 
                their positions -- seq index in seqs (chromosome order) and 
                another index at which the kmer can be found in seq
    
    The hash table data structure is applied in Python using a dictionary

    """
    hash_dict = {}
    # loop through every sequence in seqs and keep a count of iterations (i)
    for i, seq in enumerate(all_seq):
        # we use k as a step to include only NON-OVERLAPPING kmers
        for j in range(0, len(seq)-k+1, k):
            kmer = seq[j:j+k]
            try:
                # append positions to an existing kmer
                hash_dict[kmer].append((i + 1, j + 1))
            except KeyError:
                # when kmer is not yet present in kmer_dict
                hash_dict[kmer] = [(i + 1, j + 1)] 
    return hash_dict

def hitsList(query, hash_dict, k):
    """
    query : query sequence(s) inputted by the user
    hash_dict : dictionary of ref genome kmers and their positions
    k : length of kmer

    Returns
    -------
    hits_list : a sorted list containing tuples (index, shift, offset)
    
    Each sublist provides information regarding chromosome order (index),
    distance between occurance of kmer in query and ref sequence (shift) and 
    index at which the kmer can be found in ref sequence
    
    """
    hits_list = []
    # include all OVERLAPPING kmers in query
    for i in range(len(query)-k+1):
        kmer = query[i:i+k]
        # if kmer cannot be found in hash_dict, continue to NEXT iteration
        if not kmer in hash_dict:
            continue
        # extract position of kmer from hash table
        pos = hash_dict[kmer]
        for j in pos:
            # shift = kmer position in ref sequence - i (current position in query)
            shift = j[1] - i
            # index, shift, onset
            hits_list.append((j[0], shift, j[1]))  
    # sort the list by index as indicated by x[0] then by shift x[1]
    hits_list = sorted(hits_list, key = lambda x: (x[0], x[1]))
    return hits_list

#------------ALTERNATIVE CODE (IF IMPORTING LIBRARY IS ALLOWED)---------------
#def longest_match(hits_list, query_seq):
#    from itertools import groupby
#    res = [list(j) for i,j in groupby(
#    sorted(hits_list, key = lambda x: x[0] and x[1]), 
#    lambda x: x[0] and x[1])]
#    longest_match = max(res, key = len)
#    return longest_match

def masterDict(hits_list):
    """
    hits_list : a list of (index, shift, offset) of query 
                kmers that occur in hash table

    Returns
    -------
    master_dict : nested dictionary {index: {shift: [offset, ..., ...], [offset2, ...]}}
    
    """
    master_dict = {}
    for i in hits_list:
        index = i[0]
        shift = i[1]
        offset = i[2]
        # if index is not yet present in master_dict
        if not index in master_dict:
            master_dict[index] = {}   
        try:
            master_dict[index][shift].append(offset)
        except:
            # when shift is not yet present in second layer of master_dict
            master_dict[index][shift] = [offset]
    return master_dict

def longestMatch(master_dict, k):
    """
    master_dict : nested dictionary of index, shift and offset
    k : length of kmer

    Returns
    -------
    max_len : length of the longest match (int)
    long : a list of list containing index (chromosome order), 
           shift, start position and end position in ref sequence
    
    For every iteration, extract current offset and calculate current length
    by substracting the last element in offset list as indicated by 
    (offset[-1] + length of kmer) from the first element in the list.
    If current length is longer than max_len, update max_len value.
    """
    max_len = 0 
    all_length = []
    all_pos = []
    long = []
    # iterate the first layer of master_dict (index)
    for i in master_dict:
        # iterate the second layer of master_dict (shift)
        for j in master_dict[i]:
            # extract current offset
            offset = master_dict[i][j]
            curr_len = offset[-1] + k - offset[0]
            all_length.append(curr_len)
            # index, shift, start and end position
            curr_seq = [i, j, offset[0], offset[-1] + k]
            all_pos.append(curr_seq)
            if curr_len > max_len:
                max_len = curr_len
                
    for i in range(len(all_length)):
        if all_length[i] == max_len:
            long.append(all_pos[i])
            if len(long) >= 3:
                break
            
    return max_len, long
                        
def connectEnsembl(species, longest_match, Chr):
    """
    IMPORTANT: Connection cannot be established using Ubuntu 20.04 LTS,
               please use earlier version
    
    species : species name, must be specified by user
    longest_match : [index, shift, start pos, end pos]
    Chr: chromosome

    Returns gene_id and descriptin from Ensembl (e.g. AT1G03910 and 
    Cactin [Source:UniProtKB/Swiss-Prot;Acc:F4I2J8])
    
    more info: rest.ensembl.org/documentation/info/overlap_region
    """
    
    # extract start and end position of longest match in ref sequence
    start = str(longest_match[2])
    end = str(longest_match[3])
    server = "https://rest.ensembl.org"
    ext = "/overlap/region/"+species+"/"+str(Chr)+":"+start+"-"+end+"?feature=gene;" 
    # request data in json format from the server
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})  
    
    if not r.ok:
        r.raise_for_status()
        sys.exit() 
    res = r.json()
    
    gene_id = 'None'
    description = 'None'
    # if the longest match position overlaps with a gene
    if len(res) != 0:
        gene_id = res[0]["gene_id"]
        description = res[0]["description"]
    return gene_id, description

def createReport(longest_match, order, max_len, 
                  label, query, query_label, species, 
                  direction, file, screen):
    """
    longest_match : [index, shift, start, end]
    order : current index of iteration in all_query (see next function)
    max_len : length of longest match (int)
    label : a list of headers of sequences in ref genome (chromosome info)
    query : current query sequence in all_query (see next function)
    query_label : a list of headers of sequences in query file
    species : species name (str)
    direction : F (forward) or R (reverse), default value is F is no 
                specification from user
    
    Add whitespaces before and after longest match
    GCATCGACGATGCGAGC
       ||||||
       TCGACG
    """
    match = []
    longest_match_seq = []
    
    # start position in ref genome - shift = start position in query sequence
    start = longest_match[2] - longest_match[1]
    end = start + max_len
    
    # create alignment (query sequence and longest match) 
    for i in range(len(query)):
        # append whitespace 
        if i < start or i >= end:
            match.append(" ")
            longest_match_seq.append(" ")
        else:
            match.append("|")
            longest_match_seq.append(query[i])
    
    # extract the first 2 letters in label (chromosome indicator)
    Chr = label[longest_match[0]-1][:2]
    # call ensembl function to get gene id and its description 
    gene, desc = connectEnsembl(species, longest_match, Chr)
    
    # Output for two types of query: text file or direct paste to cmd
    try:
        job = query_label[order]
        # create output table
        output = tabulate([[job, str(len(query)), max_len, Chr, 
                         str(longest_match[2]), str(longest_match[3]), 
                         direction, gene, desc]], 
                         ["No", "Bases", "Longest Match", "Chr", "Start", 
                         "End", "Direction", "Gene", "Description"], 
                         tablefmt="grid")
        
    # sequence pasted directly has no label nor order
    except:
        output = tabulate([[str(len(query)), max_len, Chr, 
                         str(longest_match[2]), str(longest_match[3]), 
                         direction, gene, desc]], 
                         ["Bases", "Longest Match", "Chr", "Start", 
                         "End", "Direction", "Gene", "Description"], 
                         tablefmt="grid")
    
    # save into a text file
    with open(file+".txt", "a") as f:
        print(output, file = f)
        print(query, file = f)
        print("".join(match), file = f)
        print("".join(longest_match_seq), file = f)
        print("\n\n", file = f)
    
    # if user chooses to print the output to screen
    if screen is True:
        print(output)
        print(query)
        print("".join(match))
        print("".join(longest_match_seq))
        print("\n\n")

def main():
    """
    main() handles command line arguments and calls other functions
    
    User can input two types of query (-q or -Q)
    
    """
    #---------command line-----------
    parser = optparse.OptionParser()
    
    parser.add_option("-f", dest = "ref",
                      type = "string",
                      help = "Specify a reference genome (FASTA format).")
    
    parser.add_option("-n", dest = "name",
                      type = "string",
                      help = "Use underscore _ to separate genus and species (e.g. Arabidopsis_thaliana).")
    # query in a text file
    parser.add_option("-q", dest = "query",
                      type = "string",
                      help = "Specify a text file containing list of short sequences (e.g. query.txt).")
   # query pasted directly to command line
    parser.add_option("-Q", dest = "seq",
                      type = "string",
                      help = "Directly paste a short sequence")
    
    parser.add_option("-d", dest = "direction",
                      type = "choice",
                      choices = ["F", "R"],
                      default = "F",
                      help = "Specify orientation: F (forward) or R (reverse). The default value is F.")
    
    parser.add_option("-k", dest = "kmer",
                     type = "int",
                      default = "5",
                      help = "Specify length of kmer. The default value is 5.")
    
    parser.add_option("-o", dest = "output",
                      type = "string",
                      default = "SSAHA_report",
                      help = "Specify file name for output. The default value is SSAHA_report")
    
    parser.add_option("-s", dest = "print",
                      action = "store_true",
                      default = False,
                      help = "Option to print the result to screen")

    (options, args) = parser.parse_args() 
    
    if (options.ref is None) or (options.name is None) or ((options.query is None) and (options.seq is None)):
        parser.error("Reference genome, species name and query file must be specified")
 
    # if ref genome provided is not in fasta format
    if not any(i in options.ref for i in [".fa", ".fasta", ".fsa"]):
        parser.error("Only FASTA file (.fa, .fasta, .fsa) is accepted")
        
    ref_genome = options.ref
    ref_species = options.name
    query_file = options.query
    kmer_length = options.kmer
    orientation = options.direction
    file_name = options.output
    to_screen = options.print
 
    #--------call other functions-------------
    # parse ref genome
    ref_label, all_ref = parseFasta(ref_genome)
    hash_table = hashTable(all_ref, kmer_length)
    
    # if user input a text file query
    if not (options.query is None):
        query_label, all_query = parseFasta(query_file)
    
        for order, query_seq in enumerate(all_query):
            if orientation == "R":
                query_seq = reverseComplement(query_seq)
            hits = hitsList(query_seq, hash_table, kmer_length)
            if len(hits) != 0:
                master_list = masterDict(hits)
                max_length, longest = longestMatch(master_list, kmer_length)
                # create a report for every longest match
                for i in longest:
                    report = createReport(i, order, max_length, ref_label, query_seq, query_label, ref_species, orientation, file_name, to_screen)
                    report
            else:
                print("No hits found in database for the query")
    else: 
        query_seq = options.seq
        order = None
        query_label = None
        if orientation == "R":
            query_seq = reverseComplement(options.seq)
        hits = hitsList(query_seq, hash_table, kmer_length)
        if len(hits) != 0:
            master_list = masterDict(hits)
            max_length = longestMatch(master_list, kmer_length)
            for i in longest:
                report = createReport(i, order, max_length, ref_label, query_seq, query_label, ref_species, orientation, file_name, to_screen)
                report
        else:
            print("No hits found in database for the query")
        
if __name__ == '__main__': 
    main() 

# output memory usage and processing time
#toc = time.clock()
#print("Memory usage (in byte): "+str(process.memory_info().rss))
#print("Processing time: "+str(toc-tic))
