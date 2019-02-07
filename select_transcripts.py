#!/usr/bin/env python
#last edit: 07.02.19
#GENOMICS
import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="%(prog) reduces one FASTA file to the IDS from another file eg. to extract MAKER transcripts for proteins of a gene family from the whole genome transcript file.")

parser.add_argument("-i", dest="idfile", help="File containing sequence IDs")
parser.add_argument("-f", dest="fastafile", help="FASTA file containing sequences")
parser.add_argument("-m", dest="mode", help="mode of sequence selection. default = normal", choices=["invert", "normal"], default="normal")
parser.add_argument("-g", dest="glob", help="Specify if IDs need to be perfect matches are just have to be conatined in the sequence id. default = perfect", choices=["include", "perfect"], default="perfect")

if len(sys.argv)<2:
		parser.print_help()
		sys.exit()
Args = parser.parse_args()
if Args.fastafile == None or Args.idfile == None:
	sys.stderr.write("Error: Input file missing.")
	parser.print_help()
	sys.exit()
		
Source_name = Args.fastafile
Query_name = Args.idfile
what = Args.mode
match = Args.glob

Source = open(Source_name, "r")
Query = open(Query_name,"r")

queries = []
for line in Query:
	queries.append(line.strip())

sequences = list(SeqIO.parse(Source, "fasta"))


i=0
if what == "invert":
	sys.stderr.write("Info: Matches will be inverted.\n")
	list = list()
	if match == "perfect":
		for sequence in sequences:
			if (sequence.id not in names):
				list.append(sequence)
				i +=1
				print (">%s" % sequence.id)
				print (sequence.seq)
	if match == "include":
		for sequence in sequences:
			for name in names:
				write = True
				if (name in sequence.id):
					write = False
			if write:
				list.append(sequence)
				i +=1
				print (">%s" % sequence.id)
				print (sequence.seq)
elif what == "normal":
	sys.stderr.write("Info: Invert is set to FALSE\n")
	if match == "perfect":
		for query in queries:
			for sequence in sequences:
				if query == sequence.id:
					print (">%s" % sequence.id)
					print (sequence.seq)
					i+=1
	if match == "include":
		for query in queries:
			for sequence in sequences:
				if query in sequence.id:
					print (">%s" % sequence.id)
					print (sequence.seq)
					i+=1


string = "Wrote " + str(i) + " sequences to file\n"
sys.stderr.write(string)
Query.close()
Source.close()
