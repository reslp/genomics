#!/usr/bin/env python
# GENOMICS
# Philipp Resl
# last changed 14.07.16

import sys
import re
from Bio import SeqIO

Info = """
Identifies longest splicing isoform in Trinity assembly. 
Flags old and new help to deal with different headers in Trinity Assembly FASTA files.
old refers to the Trinity style before V2 similar to: >comp39714_c3_seq1 len=690 path=[29918167:0-689]
new is for newer Trinity style including the most recent version: >TRINITY_DN32682_c1_g1_i1 len=394 path=[372:0-393] [-1, 372, -2]

usage: get_longest_isoform.py <sequences.fasta> -old|new
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	seq_file_name = sys.argv[1]
	flag = sys.argv[2]

seq_file = open(seq_file_name,"r")

sequence_list = list()
i = 0
for sequence in SeqIO.parse(seq_file, "fasta"):
	i += 1
	#print sequence
	sequence_list.append(sequence)
	#if i > 1000:
	#	break
seq_file.close()
if "new" in flag:
	id_list = [re.search("(.+)_i\d+",sequence_list[i].id).group(1) for i in xrange(len(sequence_list))]
elif "old" in flag:
	id_list = [re.search("(.+)_seq\d+",sequence_list[i].id).group(1) for i in xrange(len(sequence_list))]
else:
	sys.stderr.write("Unknown flag specified: " + flag +"\n")
	sys.exit(0) 
#print id_list
sys.stderr.write("Total number of contigs: "+ str(len(id_list))+ "\n") 	
id_list = set(id_list)
sys.stderr.write("Total number of contigs after selecting longest splicing isoform: "+ str(len(id_list))+ "\n"	)

#sys.exit(0)
length = 0
longest = 0
i = 0
for id in id_list:
	i += 1
	for sequence in sequence_list:
		if id in sequence.id:
			length = re.search("len=(\d+)", sequence.description)
			length = int(length.group(1))
			if length > longest:
				longest = length
				name = ">" + sequence.description
				seq = sequence.seq
	sys.stderr.write("Writing sequence: " + str(i) +"\r",) 
	print name
	print seq
	length = 0
	longest = 0


	
	