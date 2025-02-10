# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:21:41 2025

@author: JiwoonPark

## Practicing coding for bioinformatics

"""

#%%
## Counting Nucleotides in a DNA sequence

DNA_string = str(input("Enter DNA string: "))

count_DNA = {}
DNA_order = ['A', 'C', 'G', 'T']

for i in DNA_string:
    count_DNA[i] = DNA_string.count(i)

print(" ".join(str(count_DNA[i]) for i in DNA_order))


#%%
## Transcribing DNA sequences into RNA sequences

DNA_string = str(input("Enter DNA string: "))

RNA_string = DNA_string.replace("T", "U")

print(RNA_string)

#%%
## Complementing DNA sequences

DNA_string = str(input("Enter DNA string: "))

complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
rev_string = reversed(DNA_string)

print("".join(complement[i] for i in rev_string))

#%%
## Calculating GC gontent in sequences

f = open("sequence.txt", 'r')

from Bio import SeqIO

seq_dict = {}
for entry in SeqIO.parse(f, "fasta"):
    seq_dict[entry.id] = str(entry.seq)

GC_content = {}
for seq_id, seq in seq_dict.items():
    A = seq.count("A")
    C = seq.count("C")
    G = seq.count("G")
    T = seq.count("T")
    GC_content[seq_id] = (C+G)/(A+C+G+T)*100

max_GC_id = max(GC_content, key=GC_content.get)
max_GC_content = GC_content[max_GC_id]

print(max_GC_id, max_GC_content)

#%%
## Counting nucleotide differences

f = "sequence.txt"
with open(f) as file:
    seq1 = file.readline().strip()
    seq2 = file.readline().strip()

hamming_distance = 0
for i in range(0, len(seq1)):
    if seq1[i] != seq2[i]: hamming_distance += 1
    
print(hamming_distance)

#%%
##  Translating RNA into Protein

f = "codon_table.txt"

RNA_codon_str = ''
with open(f) as file:
    RNA_codon_str = ",".join(entry for line in file 
                          for entry in line.replace("\t", " ").strip().split(" ")
                          if entry)

RNA_codon_list = RNA_codon_str.split(',')
RNA_codon = dict(zip(RNA_codon_list[::2], RNA_codon_list[1::2]))

RNA_str = str(input("Enter RNA string: "))

protein_str = ''
for i in range(0, len(RNA_str),3):
    if RNA_codon[RNA_str[i:i+3]] == "Stop" : break
    protein_str += RNA_codon[RNA_str[i:i+3]]
    
print(protein_str)

#%%
## Consensus and Profile

f = "sequence.txt"

import numpy as np

seq_dict = {}
for entry in SeqIO.parse(f, "fasta"):
    seq_dict[entry.id] = str(entry.seq)
strings = list(seq_dict.values())
    
profile_matrix = np.zeros((4, len(strings[0])), dtype = int)

for strand in strings:
    for i, nucleotide in enumerate(strand):
        if nucleotide == "A":
            profile_matrix[0, i] += 1
        elif nucleotide == "C":
            profile_matrix[1, i] += 1
        elif nucleotide == "G":
            profile_matrix[2, i] += 1
        elif nucleotide == "T":
            profile_matrix[3, i] += 1
conversion = {0: "A", 1: "C", 2: "G", 3: "T"}
consensus_string = "".join([conversion[np.argmax(profile_matrix[:, i])] for i in range(len(strings[0]))])

output_file = "output.txt"
with open(output_file, "w") as f:
    print(consensus_string, file=f)
    for i, nucleotide in enumerate(["A", "C", "G", "T"]):
        print(f"{nucleotide}: {' '.join(map(str, profile_matrix[i]))}", file=f)
    
#%%
##  Finding a motif in two sequences

st1, st2 = [str(x) for x in input("Enter DNA strings 1 and 2: ").split()]

location = []
for i in range(len(st1)):
    if st1[i:i+len(st2)] == st2 : location.append(i+1)

print(' '.join(map(str, location)))

#%%
## Finding a shared motif between multiple seqeunces

from Bio import SeqIO

f = "sequence.txt"
seq_dict = {}
for entry in SeqIO.parse(f, "fasta"):
    seq_dict[entry.id] = str(entry.seq)
strings = list(seq_dict.values())

def longest_common_substring(strings):
    ref_seq = min(strings, key=len)
    seq_len = len(ref_seq)

    for n in range(seq_len, 1, -1): 
        substrings = {ref_seq[i:i+n] for i in range(seq_len - n + 1)} 
    
        for substr in substrings:
            if all(substr in seq for seq in strings):  
                return(substr) 
    return("No common substring")

print(longest_common_substring(strings))

#%%
##  Finding ORFs (Open Reading Frames)

f = "C:/Users/JiwoonPark/Downloads/RNA codon table.txt"

RNA_codon_str = ''
with open(f) as file:
    RNA_codon_str = ",".join(entry for line in file 
                          for entry in line.replace("\t", " ").strip().split(" ")
                          if entry)
RNA_codon_list = RNA_codon_str.split(',')
RNA_codon = dict(zip(RNA_codon_list[::2], RNA_codon_list[1::2]))

DNA_string = str(input("Enter DNA string: "))

complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
DNA_string_comp = ("".join(complement[i] for i in DNA_string))
DNA_string_comp_rev = DNA_string_comp[::-1]

RNA_string = DNA_string.replace("T", "U")
RNA_string_comp_rev = DNA_string_comp_rev.replace("T", "U")

def translate_rna(RNA_str):
    protein_str = []
    
    for shift in range(3):  
        i = shift

        while i < len(RNA_str) - 2:
            if RNA_str[i:i+3] == "AUG": 
                protein_seq = ""
                stop_codon = False
                for j in range(i, len(RNA_str) - 2, 3):
                    codon = RNA_str[j:j+3]
                    amino_acid = RNA_codon.get(codon, None)
                    if amino_acid == "Stop": 
                        stop_codon = True
                        break
                    if amino_acid is None: break
                    protein_seq += amino_acid        
                if stop_codon and protein_seq:
                    protein_str.append(protein_seq)
            i += 3
            
    return(protein_str)

all_strings = translate_rna(RNA_string) + translate_rna(RNA_string_comp_rev)
print("\n".join(set(all_strings)))


#%%
## RNA splicing and translation

from Bio import SeqIO

f = "sequences.txt"
seq_dict = {}
for entry in SeqIO.parse(f, "fasta"):
    seq_dict[entry.id] = str(entry.seq)
strings = list(seq_dict.values())

DNA_string, introns = strings[0], strings[1:]

for i in range(len(introns)):
    DNA_string = DNA_string.replace(introns[i], "")

RNA_string = DNA_string.replace("T", "U")

protein_str = ''
for i in range(0, len(RNA_string) - (len(RNA_string)%3),3):
    if RNA_codon[RNA_string[i:i+3]] == "Stop" : break
    protein_str += RNA_codon[RNA_string[i:i+3]]
    
print(protein_str)


#%%
## Finding a shared spliced motif 

seq1, seq2 = [str(x) for x in input("Enter DNA strings 1 and 2: ").split()]

def longest_motif(x, y):
    seq_matrix = [["" for _ in range(len(y)+1)] for _ in range(len(x) + 1)]
    
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            if x[i-1] == y[j-1]: 
                seq_matrix[i][j] = seq_matrix[i-1][j-1] + x[i-1]
            else: 
                seq_matrix[i][j] = max(seq_matrix[i-1][j], seq_matrix[i][j-1], key = len)
                
    return seq_matrix[len(x)][len(y)]

print(longest_motif(seq1, seq2))


#%%
## Identifying reverse palindromes

DNA_string = str(input("Enter DNA string: "))

complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def rev_comp(DNA_str):
    return("".join(complement[i] for i in reversed(DNA_str)))

DNA_substring, DNA_substring_rev_comp, pos_length = [], [], []

for j in (range(4,13,2)):
    for i in range(len(DNA_string)-(j-1)):
        DNA_substring.append(DNA_string[i:i+j])
        DNA_substring_rev_comp.append(rev_comp(DNA_substring[-1]))
        if DNA_substring[-1] == DNA_substring_rev_comp[-1]: 
            pos_length.append((i+1, j))
    
print("\n".join(f"{x} {y}" for x, y in pos_length))


#%%
# Calculating transition/transversion ratio

seq1, seq2 = [str(x) for x in input("Enter DNA strings 1 and 2: ").split()]
transit_count = 0
transver_count = 0

transit_dict = [("A","G"), ("G","A"), ("T","C"), ("C","T")]
transver_dict = [("A","C"), ("C","A"), ("A","T"), ("T","A"), ("T","G"), ("G","T"), ("G","C"), ("C","G")]

for i in range(0, len(seq1)):
    if len(seq1) != len(seq2): break
    for seq in zip(seq1, seq2):
        if seq in transit_dict: transit_count += 1
        if seq in transver_dict: transver_count += 1
        
print(transit_count/transver_count)
