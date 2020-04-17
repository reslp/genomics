#!/usr/bin/env python

import argparse
import sys
from collections import defaultdict


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="filter_blast_results.py", description = """This script filters BLAST results in outfmt 6 qseqid sseqid pident length slen qlen evalue bitscore field order created for the Pmem temp experiment.\n Filtering is based on the following criteria: pvalue < 1e-10""", epilog = """written by Philipp Resl""")
pars.add_argument('-i', dest="i", required=True, help="Path to the blast results file.")
args=pars.parse_args()

print("Will filter according to the following criteria:", file=sys.stderr)
print("\t1. evalue < 1e-10", file=sys.stderr)
print("\t2. coverage (length / qlen) >= 0.7 (70%)", file=sys.stderr)
print("\t3. Hit is unqiue for query (not found in hits for other queries)", file=sys.stderr)
print("\nImportant: The script expects the blast results to be formatted in non-standard tabular format:", file=sys.stderr)
print("qseqid sseqid pident length slen qlen evalue bitscore\n", file=sys.stderr)


file = open(args.i, "r")
gene_dict = defaultdict(list)
for line in file:
	line = line.strip()
	qseqid, sseqid, pident, length, slen, qlen, evalue, bitscore = line.split("\t") 
	qseqid = qseqid.split("|")[-1].split("_")[0]
	if float(evalue) <= 1e-10: # first discard all results with an e-value > 1e-10
		if float(length) / float(qlen) >= 0.7: # next only keep seqs with coverage >= 70% of the query.
			if sseqid not in gene_dict[qseqid]: # avoid double entries
				gene_dict[qseqid].append(sseqid)
				#print(qseqid, sseqid, evalue)
file.close()
# some hits could be hits fro multiple queries. It is necessary to filter them out when the occur in more than query.
protein_ids = []
for gene in gene_dict.keys():
	 for id in gene_dict[gene]:
	 		protein_ids.append(id)
	 		ids_to_exclude = []
for id in protein_ids:
	if protein_ids.count(id) > 1:
		ids_to_exclude.append(id)
ids_to_exclude = set(ids_to_exclude)
for id in ids_to_exclude:
	occurence = 0
	for gene in gene_dict.keys():
		if id in gene_dict[gene]:
			occurence += 1
	if occurence > 1:
		print(id, "was found",occurence,"times, will drop it\t",file=sys.stderr)
		for gene in gene_dict.keys():
			if id in gene_dict[gene]:
				gene_dict[gene].remove(id)

# print final selection of genes
for gene in gene_dict.keys():
	if len(gene_dict[gene]) > 0:
		for id in gene_dict[gene]:
			print("potential_", gene, "_",gene_dict[gene].index(id)+1, "\t", id, sep="")
