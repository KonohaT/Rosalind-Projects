#My solutions to the problems in https://rosalind.info/problems/list-view/

codon_table = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V", "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V", "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
"UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V", "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A", "UCC" : "S", "CCC" : "P", "ACC" : "T",
"GCC" : "A", "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A", "UCG" : "S",
"CCG" : "P", "ACG" : "T", "GCG" : "A", "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D", "UAC" : "Y", "CGU" : "R", "AGU" : "S",  "GGU" : "G", 
"CAC" : "H", "AAC" : "N", "GAC" : "D", "UAA" : "Stop", "CAA" : "Q", "AAA" : "K", "GAA" : "E", "UAG" : "Stop", "CAG" : "Q", "AAG" : "K", "GAG" : "E", "UGU" : "C",
"UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G", "UGA" : "Stop", "CGA" : "R", "AGA" : "R", "GGA" : "G", "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"}

#Problem 1
def countBases(strand):
    returnString = ""
    nucleotides = ['A', 'C', 'G', 'T']
    for i in range(len(nucleotides)):
        counted = str(strand.count(nucleotides[i]))
        returnString += counted + " "
    returnString = returnString.strip()
    return returnString

#Problem 2
def transcribeRNA(strand):
    return strand.replace("T", "U")

#Problem 3
def complementDNA(strand):
    reverseStrand = strand [::-1]

    replacementTable = {"A": "T", "C":"G", "G":"C", "T":"A"}

    finalStrand = ""
    for i in reverseStrand:
        try:
            finalStrand += replacementTable[i]
        except KeyError:
            print("One of these characters was not a valid base")
            break
    return finalStrand


#Problem 4 Fibbonacci 1
def count_Rabbits(months, offspringsPer = 1):
    months = int(months)
    offspringsPer = int(offspringsPer)

    minus2 = 0
    minus1 = 1
    current = 1
    monthsPassed = 1

    while monthsPassed < months - 1: #-1 because the monthsPassed starts at 1 in the problem
        minus2 = minus1
        minus1 = current
        current = minus1 + (minus2 * offspringsPer)

        monthsPassed += 1
    
    return current

#Problem 4.1 - Fibonacci 1 redone in NumPy
import numpy as np
def count_array_rabbits(months, offsprings_per = 1):
    totals = np.array([1,1])
    i = 2
    while i < months:
        add = totals[-1] + offsprings_per * totals[-2]
        totals = np.append(totals, add)
        i += 1
    return totals[-1]

#Problem 5 GC Content
from operator import attrgetter
def findGC(totalString): #Should be able to take entire input as one string

    #splits the string into multiple DNA strands with names
    unprocessedList = totalString.split(">")
    del unprocessedList[0]

    #defines the FASTA file format as a class
    class FASTA:
        def __init__(self, name, strand):
            self.name = name
            self.strand = strand
            self.gcContent = (strand.count("G") + strand.count("C")) / len(strand) * 100 #percentage
        #defines the gcCount as the number of Gs and Cs

    #instantiates each FASTA instance, separating name and strand.
    processedList = []
    for instance in unprocessedList:
         split = instance.split("\n")
         processedList.append(FASTA(split[0], ''.join(split[1:])))

    #compares gcContent to find the largest and returns that one. 
    maxName = max(processedList, key = attrgetter('gcContent')).name
    maxGCContent = max(processedList, key = attrgetter('gcContent')).gcContent    
    print(maxName)
    print(maxGCContent)

#Problem 6 - Count Point Mutations
def countPointMuts(strand1, strand2):
    muts = 0

    if len(strand1) != len(strand2):
        print("Strands are of different lengths!")
        return 1
    
    for basePos in range(len(strand1)):
        if strand1[basePos] != strand2[basePos]:
            muts += 1
    
    return muts

#Problem 7
def find_mendelian(hom_dom, hetero, hom_rec):
    total_pop = hom_dom + hetero + hom_rec


#Problem 8 - mRNA to Amino Acids
def translate_RNA(strand):
    protein_strand = ""

    starting_index = strand.find("AUG")
    if starting_index == -1:
        raise Exception("No start codon found")
    
    i = starting_index
    while i < len(strand):
        current_codon = strand[i:i+3]      
        if codon_table[current_codon] == "Stop":
            return protein_strand
        else:
            try:
                protein_strand += codon_table[current_codon]
            
            except KeyError:
                raise Exception("Not an acceptable codon")
        i += 3

    return protein_strand
#Problem 9 Finding a motif
def find_motif(strand, sub_strand):
    start_pos = 0
    found_pos = 0
    found = []
    while True:
        found_pos = strand.find(sub_strand, start_pos, len(strand))
        start_pos = found_pos + 1

        if start_pos > len(strand) or found_pos == -1:
            break
        found.append(found_pos + 1)

    for i in found:
        print(i, end = " ")
    print("")

#Problem 10 Finding Mortal Rabbits
def count_mortal_rabbits(month, lifespan):
    totals = [1,1]
    pos = 2
    while pos < month:
        if pos < lifespan:
            totals.append(totals[-1] + totals[-2])
        elif pos == lifespan:
            totals.append(totals[-1] + totals[-2] - 1)
        else:
            totals.append(totals[-1] + totals[-2] - totals[-lifespan-1])
        pos += 1
    return totals[-1]

#Problem 11 - Expected Offsprings
from numpy import sum, array
from math import floor
def calc_dom_phen(pairings: list):
    if not isinstance(pairings, list):
        raise TypeError
    offsprings_per = 2
    arr = array(pairings)
    probabilities = array([1,1,1,0.75, 0.5, 0])
    return(sum(arr * probabilities) * offsprings_per)

from Bio import SeqIO
import numpy as np
def find_consensus(filename):
    if not isinstance(filename, str):
        raise TypeError
    
    sequences = []

    handle = open(filename, 'r')

    
    for record in SeqIO.parse(handle, 'fasta'):
        sequence = []
        for nt in record.seq:
            sequence.extend(nt)
        sequences.append(sequence)
    
    strand_array = np.array(sequences)
    nt_list = ["A", "C", "G", "T"]
    profile = np.zeros([4, len(sequences[0])], dtype=str) #Profile order is A,C,G,T
    consensus = np.zeros(len(sequences[0]), dtype=str)

    for column in range(len(strand_array[0])):
        i = 0
        for nt in nt_list:
            profile[i, column] = str(np.count_nonzero(strand_array[:, column] == nt))
            i += 1

        elements, counts = np.unique(strand_array[:, column], return_counts=True)
        consensus[column] = elements[np.argmax(counts)]

    print("".join(consensus))

    for i in range(len(profile)):
        print(nt_list[i] + ": " + " ".join(profile[i]))

from Bio import SeqIO
import requests
import re
from time import sleep

def find_nglyco_pos(cID):
    if not isinstance(cID, str):
        raise TypeError
    clean_cID = cID.split("_")[0]

    #get data from url
    url = "http://www.uniprot.org/uniprot/" + clean_cID + ".fasta"
    fasta_output = requests.get(url).text

    #FASTA processing
    seq = fasta_output.split('\n', 1)[1].replace('\n', '')
    
    #Regex search for positions
    pattern = r'(?=(N[^P][ST][^P]))'
    matches = []
    for match in re.finditer(pattern, seq):
        matches.append(str(match.start()+1))
    
    if matches == []:
        return 0

    #Clean output formatting
    print(cID)
    print(" ".join(matches))
    sleep(1)

def find_all_nglyco(id_list):
    for id in id_list:
        find_nglyco_pos(id)

def infer_mRNA(prot_strand):
    possibilities_table = {"A": 4, "C":2, "D":2, "E":2, "F":2, "G":4, "H":2, "I":3, "K":2, "L":6, "M":1, "N":2, "P":4, "Q":2, "R":6, "S":6, "T":4, "V":4, "W":1, "Y":2}
    
    permutations = 1
    for prot in prot_strand:    
        if prot not in possibilities_table:
            return TypeError
        else:
            permutations *= possibilities_table[prot]
    permutations *= 3


    return permutations % 1000000

"""
import re
from Bio import SeqIO
def find_protein_from_start(seq, start_location):
    prot_seq = ""
    current_pos = start_location

    while True:
        
        try:
            codon = seq[current_pos:current_pos+3]
        except: //fix this section
        
        if codon_table[codon] == "Stop":
            break
        elif codon not in codon_table:
            print("Reached unreadable codon")
            break
        else:
            prot_seq.append(codon_table[codon])
            current_pos += 3

def find_ORFs(filename): #Data in the .fasta file
    

    handle = open(filename, 'r')

    for record in SeqIO.parse(handle, 'fasta'):
        seq = str(record.seq)
    print(seq)

    start = r'(?=ATG)'
    start_locations = list(match.start() for match in re.finditer(start, seq))
    print(start_locations)

    for location in start_locations:
            


def find_permut(highest_int):
    if not isinstance(highest_int, int):
        raise TypeError

find_ORFs("sampledata.fasta")
"""


def find_ORF2(): #regex not working, debug later...
    from Bio import SeqIO
    from pathlib import Path
    import re
    data_path = Path("Rosalind-Projects/data/")
    file_name = data_path / "rosalind_orf.txt"
    file = open(file_name, 'r')
    for segment in SeqIO.parse(file, 'fasta'):
        sequence = str(segment.seq)
    
    #ORF_pattern = re.compile("A[TU]G(...)*(?:[TU]AA|[TU]AG|[TU]GA)")
    ORF_pattern = re.compile('(?=(ATG(?:...)*?)(?=TAA|TAG|TGA))')
    ORFs = re.match(ORF_pattern, sequence)
    print(ORFs)

    """
    def find_start_locations(sequence):
        start_codon = re.compile("A[TU]G")
        start_locations = list(match.start for match in re.finditer(start_codon, sequence))
        return start_locations
    
    print(find_start_locations(sequence))
    """

def finding_a_shared_motif_alt():
    #This would solve all of my problems, but it's no fun: https://pypi.org/project/pylcs/
    from Bio import SeqIO
    from pathlib import Path
    from itertools import combinations
    data_path = Path("Rosalind-Projects/data/")
    file_name = data_path / "rosalind_orf.txt" #replace
    seq_list = []
    file = open(file_name, 'r')
    for segment in SeqIO.parse(file, 'fasta'):
        sequence = str(segment.seq)
        seq_list.append(sequence)
    
    root = {"base": "R", "children": []}

    def find_all_substrings(string):
        length = len(string) + 1
        return [string[x:y] for x, y in combinations(range(length), r=2)]

    def create_child_node(parent, base):
        new = {"base": base, "children": []}
        parent["children"].append(new)
        return new

    def add_new_substring(substring):
        current = root
        for known_base in substring:
            child_bases = [x["base"] for x in current["children"]]
            if known_base not in child_bases:
                current = create_child_node(current, known_base)
            else: 
                current = next(x for x in current["children"] if x["base"] is known_base)

    #def find_longest_branch(root):


    example_string = "ATGGC"
    substrings = find_all_substrings(example_string)
    for substring in substrings:
        add_new_substring(substring)
    return root


def finding_a_shared_motif():
    from Bio import SeqIO
    from pathlib import Path
    from itertools import combinations
    data_path = Path("Rosalind-Projects/data/")
    file_name = data_path / "rosalind_lcsm.txt"
    seq_list = []
    file = open(file_name, 'r')
    for segment in SeqIO.parse(file, 'fasta'):
        sequence = str(segment.seq)
        seq_list.append(sequence)
    
    def find_all_substrings(string):
        length = len(string) + 1
        return [string[x:y] for x, y in combinations(range(length), r=2)]
    
    substrings_of_first = find_all_substrings(seq_list[0])
    common_substrings = []
    for substring in substrings_of_first:
        if all(substring in string for string in seq_list[1:]):
            common_substrings.append(substring)

    return max(common_substrings, key=len)

def genbank_intro():
    from pathlib import Path
    from Bio import Entrez
    data_path = Path("Rosalind-Projects/data/")
    file_name = data_path / "rosalind_gbk.txt"
    file = open(file_name, 'r')
    genus_name = file.readline().strip()
    date_1 = file.readline().strip()
    date_2 = file.readline().strip()

    Entrez.email = "konohatomonoduval@gmail.com"
    #handle = Entrez.esearch(db="nucleotide", term=genus_name + "[Organism]")
    handle = Entrez.esearch(db="nucleotide", term=genus_name + "[Organism]", mindate=date_1, maxdate=date_2, datetype="pdat")
    record = Entrez.read(handle)
    print(date_1 + " - " + date_2)
    return record["Count"]

def data_formats():
    from pathlib import Path
    from Bio import Entrez, SeqIO
    data_path = Path("Rosalind-Projects/data/")
    file_name = data_path / "rosalind_frmt.txt"
    file = open(file_name, 'r')
    ids = file.readline().strip()
    ids = ids.split(" ")

    Entrez.email = "konohatomonoduval@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta")
    record = list(SeqIO.parse(handle, 'fasta'))
    
    current_length = len(record[0].seq)
    for segment in record:
        if len(segment.seq) < current_length:
            current_length = len(segment.seq)
            shortest = segment

    return shortest.format("fasta")

