#!/usr/bin/env python
# written by Philipp Resl
#GENOMICS
import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Outputs part of a sequence based on start and end position for a given sequence ID of a sequence file in FASTA format.")

parser.add_argument("-i", "--id", dest="id", help="Sequence ID.")
parser.add_argument("-f", "--fasta", dest="fastafile", help="FASTA file containing sequences.")
parser.add_argument("-s", "--start", dest="start", action="store", help="Start Position in sequence.", type = int, default=0)
parser.add_argument("-e", "--end", dest="end", action="store", help="End position in sequence.", type = int, default =0)
Args = parser.parse_args()

if len(sys.argv)<2:
	parser.print_help()
	sys.exit()
print(Args)
if Args.fastafile == None:
	sys.stderr.write("Error: Input file missing.")
	parser.print_help()
	sys.exit()
	
Source_name = Args.fastafile

Source = open(Source_name, "r")
sequences = list(SeqIO.parse(Source, "fasta"))
Source.close()

i=0
for sequence in sequences:
	if (Args.id == sequence.id.split(">")[-1]):
		i +=1
		if Args.start > Args.end: #inverted hits
			start = Args.end
			end = Args.start
		else:
			start = Args.start
			end = Args.end
		
		if start - 1 < 0:
			print("Start is <0. Abort")
			sys.exit(1)
		if end + 1 > len(sequence.seq):
			print("End position larger than sequence length. Abort")
			sys.exit(1)

		print (">%s_start:%s_end:%s" % (sequence.id, start, end))
		print (sequence.seq[start-1:end])

