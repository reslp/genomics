#!/usr/bin/env python
#last edit: 29.06.16
#GENOMICS
import sys
from Bio import SeqIO
import pandas as pd

info = """
This script counts occurences of functional annotation terms for a blast2go output file for sequences a given fastafile.

Usage:
blast2go_stat.py <bast2gofile> <sequencefile> [-go | -ip]

<infile> is a blast2go tab formated output file with columns named like:
SeqName	Description	Length	#Hits	e-Value	sim_mean	#GO	GO_Names_list	Enzyme_Codes_list	InterPro_IDs

flags [-go | -ip] lets you extract either GO terms or InterPro IDs.

Output may be piped into an output file.

"""

if len(sys.argv) < 3:
	sys.stderr.write(info)
	quit()
else:
	file_name = sys.argv[1]
	seq_file_name= sys.argv[2]
	what = sys.argv[3]

sys.stderr.write("Reading sequence file... \n")
seq_file = open(seq_file_name,"r")
seq_list = list()
for sequence in SeqIO.parse(seq_file,"fasta"):
	seq_list.append(sequence.id)
seqlistlen = str(len(seq_list))
sys.stderr.write("Found " + seqlistlen + " sequences \n")
seq_file.close()

#file = open(file_name, "r")
DB = pd.read_csv(file_name,sep="\t")
DB = DB.set_index("SeqName")
#print DB.index.values.tolist()
terms = list()
length = len(DB.index.values.tolist())
fin_len = 0
if what == "-go":
	sys.stderr.write("Extracting GO terms...\n")
	for sequence in seq_list:
		#print sequence
		if sequence in DB.index.values.tolist():
			fin_len += 1
			terms.append(str(DB.ix[sequence].GO_Names_list).split("; "))
	terms = [val for sublist in terms for val in sublist] #flatten list
	terms = [item for item in terms if item != "nan"] #remove missing values
	terms_list = list()
	for term in set(terms):
		terms_list.append((term, terms.count(term)))
	terms_list = sorted(terms_list, key=lambda tup: tup[1]) #sorted
	for item in terms_list:
		print item[0], "\t", item[1], "\t", float(item[1])/fin_len
	sys.stderr.write("Total number of proteins in Blast2Go file: " + str(length) +"\n")
	sys.stderr.write("Number of proteins in selected set: " + seqlistlen +"\n")
	sys.stderr.write("Found GO terms for number of sequences: " + str(fin_len) +"\n")
	sys.stderr.write("Total number of GO Terms: " + str(len(terms_list)) +"\n")
elif what == "-ip":
	sys.stderr.write("Extracting InterPro IDs...\n")
	for sequence in seq_list:
		if sequence in DB.index.values.tolist():
			fin_len += 1
			terms.append(str(DB.ix[sequence].InterPro_IDs).split("; "))
	terms = [val for sublist in terms for val in sublist] #flatten list
	terms = [item for item in terms if item != "no IPS match"]
	terms = [item.split(" ")[0] for item in terms]
	#print terms
	#quit()
	terms_list = list()
	for term in set(terms):
		terms_list.append((term, terms.count(term)))
	terms_list = sorted(terms_list, key=lambda tup: tup[1]) #sorted
	for item in terms_list:
		print item[0], "\t", item[1], "\t", float(item[1])/fin_len
	sys.stderr.write("Total number of proteins in Blast2Go file: " + str(length) +"\n")
	sys.stderr.write("Number of proteins in selected set: " + seqlistlen +"\n")
	sys.stderr.write("Found InterPro terms for number of sequences: " + str(fin_len) +"\n")
	sys.stderr.write("Total number of Interpro IDs: " + str(len(terms_list)) +"\n")

else:
	sys.stderr.write("Please specify what to extract.\n")
