#!/usr/bin/env python
#last edit: 04.09.2016

import sys
import operator
#import pandas as pd #needs >0.16
from Bio import SeqIO

info = """
Filters BLAST tab output for sequences without hits in other species according to evalue cutoff.

Usage:
get_unique.py <blast_tab_inputfile> <sequence_file.fas>

"""

if len(sys.argv) < 2:
	sys.stderr.write(info)
	quit()
else:
	blast_file_name = sys.argv[1]
	sequence_file_name  = sys.argv[2]

evalue = 1e-2
#alen = 0.5 #= 50%


#setdefault(t,[])
blast_file = open(blast_file_name, "r")
evalue_dict= {}
names_list = []
for line in blast_file:
	names_list.append(line.split("\t")[0])
	evalue_dict.setdefault(line.split("\t")[0],[]).append(line.split("\t")[10])
blast_file.close()

names_list = set(names_list)
unique_list = []
for name in names_list:
	if evalue < float(min(evalue_dict[name], key=operator.itemgetter(1))):
		unique_list.append(name)
		continue

print "There is a total of",str(len(unique_list)),"unique genes without blast hits with an evalue under", str(evalue)


print "Reading sequence file"
seq_file = open(sequence_file_name,"r")
sequences = SeqIO.parse(seq_file, "fasta")
sequence_list = [sequence for sequence in sequences] #memory inefficient!
seq_file.close()

print "Writing output files:"
species_list = ["agyrium", "lambiella", "trapelia", "xylographa", "graphis", "rimularia", "trapeliopsis", "icmadophila", "dibaeis"]
for species in species_list:
	file = open(species+"_unqiue_genes.fas","w")
	i = 0
	for sequence in sequence_list:
		if species in sequence.id and sequence.id in unique_list:
			i += 1
			file.write(">"+sequence.id+"\n")
			file.write(str(sequence.seq)+"\n")
			#sequence_list.
	print "Wrote", str(i), "sequences for", species
	file.close()
