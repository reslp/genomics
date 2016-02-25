#!/usr/bin/env python
#last edit: 24.02.16
#GENOMICS
import sys
from Bio import SeqIO

Info = """Reduced one FASTA file to the IDS from another eg. to extract MAKER transcripts for proteins of a gene family from the whole genome transcript file.

usage: select_transcripts.py <sourcefile> <queryfile> <invert|normal> > outfile

<invert|normal> decides if script should invert selection or not.

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
Query = open(Query_name,"r")

queries = list(SeqIO.parse(Query, "fasta"))
sequences = list(SeqIO.parse(Source, "fasta"))

i=0
if what == "invert":
	sys.stderr.write("Info: Invert is set to TRUE\n")
	names = list()
	list = list()
	for query in queries:
		names.append(query.id)
	for sequence in sequences:
		if (sequence.id not in names):
			list.append(sequence)
			i +=1
			print ">%s" % sequence.id
			print sequence.seq
elif what == "normal":
	sys.stderr.write("Info: Invert is set to FALSE\n")
	for query in queries:
		for sequence in sequences:
			if query.id == sequence.id:
				print ">%s" % sequence.id
				print sequence.seq
				i+=1
else:
	sys.stderr.write("Please specify either normal or invert\n")
	sys.stderr.write(Info)
	quit()


string = "Wrote " + str(i) + " sequences to file\n"
sys.stderr.write(string)
Query.close()
Source.close()
