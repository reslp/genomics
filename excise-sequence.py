#!/usr/bin/env python3
import argparse


def extract_sequences(infile, seqid, positions):
	with open(infile, "r") as inf:
		found = False
		seq = ""
		for line in inf:
			if line.startswith(">") and line.strip() == ">"+seqid:
				found = True
				continue
			if found and not line.startswith(">"):
				seq += line.strip()
	for pos in positions:
		start, stop = pos.split(",") 
		print(">"+seqid+"_"+start+"_"+stop)
		print(seq[int(start)-1: int(stop)])


if __name__ == '__main__':
	pars = argparse.ArgumentParser(prog="excise-sequence.py", description = """This script excise a part of a longer sequence from a FASTA file.""", epilog = """written by Philipp Resl""")
	pars.add_argument('--in', dest="inf", action="store", help="Path to the input file in FASTA format.", required=True)
	pars.add_argument('--id', dest="id", action="store", help="The sequence ID as it is in the FASTA file (without the >)", required=True)
	pars.add_argument('--pos', dest="pos", action="append", help="The sequence to be excised defined by start and stop position seperated by comma. Multiple arguments are possible. Eg. --pos 1,25 --pos 267,5340", required=True)
	args=pars.parse_args()

	extract_sequences(args.inf, args.id, args.pos)
		


