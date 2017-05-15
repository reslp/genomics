#!/usr/bin/env python
#last edit: 29.06.16
#GENOMICS
import sys
#import operator
from collections import defaultdict

info = """
This script filters PFAM domains from files created with pfam.xfam.org

Usage:
filter_pfam.py <inputfile>

"""

if len(sys.argv) < 2:
	sys.stderr.write(info)
	quit()
else:
	file_name = sys.argv[1]


file = open(file_name, "r")

pfam_dict = defaultdict(list)
for line in file: 
	gene_name = line.split("\t")[0]
	pfam_name = line.split("\t")[5]
	pfam_dict[pfam_name].append(gene_name)
	
pfam_tuples = [(key, pfam_dict[key]) for key in pfam_dict.keys()]
pfam_sorted = sorted(pfam_tuples, key = lambda x: len(x[1]))

for term in pfam_sorted:
	gene_string = ""
	for item in term[1]:
		gene_string += item 
		gene_string += " "
	print term[0] + "\t" + str(len(term[1])) +"\t"+ gene_string
