#!/usr/bin/env python3
# this has to be run with uv run scriptname for the dependencies to be installed and used!
# /// script
# dependencies = [
#   "biopython",
# ]
# ///


import sys

import argparse
import sys
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="filter_blast_results.py", description = """Filter blast and extract sequences from hits""", epilog = """written by Philipp Resl""")
pars.add_argument('-b', dest="b", required=True, help="Path to the blast results file.")
pars.add_argument("-g", dest="g", required=True, help="Path to genome file for extraction.")
pars.add_argument("--id", dest="id", help="Gene ID to append to sequence name.", default="")
args=pars.parse_args()

blast_hits = []
with open(args.b, "r") as blastfile:
	for line in blastfile:
		if "\t0.0\t" in line: # only keep hits with e value of 0.0
			blast_hits.append(line.strip())

contigs_with_hits = []
for hit in blast_hits:
	contigs_with_hits.append(hit.split("\t")[1])
contigs_with_hits = set(contigs_with_hits)

hit_ranges = {}
for contig in contigs_with_hits:
	#print("Hits for", contig)
	start = []
	for hit in blast_hits:
		if contig in hit:
			#print(hit)
			start.append(int(hit.split("\t")[8]))
			start.append(int(hit.split("\t")[9]))
	#start = start.sort()
	start = sorted(start)
	hit_ranges[contig] = [start[0], start[-1]] # creates a 50bp overhang on both sides...


for record in SeqIO.parse(args.g, "fasta"):  
	for contig in hit_ranges.keys():
		if contig == record.id:
			print(">"  + record.id +"_hit_" + str(hit_ranges[contig][0]) + "-" + str(hit_ranges[contig][1])+"_"+args.id)
			print(record.seq[hit_ranges[contig][0]: hit_ranges[contig][1]])







