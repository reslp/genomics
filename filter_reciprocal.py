#!/usr/bin/env python
# written by Philipp Resl
# 24.08.16

import sys
#from Bio import SearchIO
from collections import defaultdict

info = """
Filters BLAST results for reciprocal best hits for a group of species.
Input file needs to be tab-delimited blast file reduced to only best hits for each species with get_best_hit.py

Usage:
filter_reciprocal.py <blast_results.txt> <nspecies>

"""

if len(sys.argv) < 2:
	sys.stderr.write(info)
	quit()
else:
	file_name = sys.argv[1]
	nspecies = int(sys.argv[2])

file = open(file_name,"r")
genes = defaultdict(list)
for line in file:
	genes[line.split("\t")[0]].append(line.split("\t")[1])
file.seek(0)

redundant = []
keys = genes.keys()
keys = sorted(keys)

reciprocs = []
j= 0
for gene in keys:
	i = 0
	output = ""
	#print genes[gene]
	for hit in sorted(genes[gene]):
		#print hit
		if sorted(genes[gene]) == sorted(genes[hit]):
			i += 1
			#output +=  hit + "\t"
		if i == nspecies:
			j += 1
			#print gene
			reciprocs.append(sorted(genes[gene]))
			break
#print j	
#print len(reciprocs)
final_list = []
for sublist in reciprocs:
	if sublist not in final_list:
		final_list.append(sublist)
#print len(final_list)
output = ""
for list in final_list:
	for element in list:
		output += element + "\t"
	print output
	output = ""
		

#print len(final_list)
	
	
	
	