#!/usr/bin/env python
#last edit: 15.06.16
#GENOMICS


import sys
from Bio import SeqIO
import pandas as pd

Info = """Reduced a FASTA file to the IDS from dbCAN output according to a specific group of cazymes"

usage: select_cazys.py <sourcefile.fas> <queryfile.txt> CAZYGROUP > outfile

CAZYGROUP is a CAZy group ID eg. GH, AA, CE
"""
invert = True
if len(sys.argv) < 4:
	sys.stderr.write(Info)
	quit()
else:
	Source_name = sys.argv[1]
 	Query_name = sys.argv[2]
 	what = sys.argv[3]

Source = open(Source_name, "r")

queries = pd.read_csv(Query_name,sep="\t", header=None)
sequences = list(SeqIO.parse(Source, "fasta"))

for index, query in queries.iterrows():
	for sequence in sequences:
		if (query[0] in sequence.id) and (what in query[1]):
			print ">"+ sequence.id
			print sequence.seq
			




