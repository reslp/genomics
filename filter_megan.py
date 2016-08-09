#!/usr/bin/env python
# GENOMICS
# Philipp Resl
# last changed 08.07.16

import sys
from Bio import SeqIO

Info = """
Filters FASTA files according to a tab delimited file from MEGAN (export CSV-> readName_to_taxonPath).

usage: filter_megan.py <sequences.fasta> <megan_file> <<taxonomic_level>>
"""

if len(sys.argv) < 4:
	sys.stderr.write(Info)
	quit()
else:
	seq_file_name = sys.argv[1]
	megan_file_name = sys.argv[2]
	tax_level = sys.argv[3]

seq_file = open(seq_file_name,"r")
megan_file = open(megan_file_name,"r")

records = list()
for record in megan_file:
	records.append(record)
megan_file.close()

for sequence in SeqIO.parse(seq_file, "fasta"):
	for record in records:
		record_split = record.split("\t")
		if tax_level in record_split[1] and sequence.id == record_split[0]:
			print ">" + sequence.id
			print sequence.seq	
			break		
seq_file.close()

		
