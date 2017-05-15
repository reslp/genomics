#!/usr/bin/env python
# GENOMICS
# Philipp Resl
# last changed 12.01.17

import sys
import time
Info = """
Replaces GO Term Names with GO ID Numbers. Creates file that will work with the topGO package as a gene-to-GO file.

Takes a BLAST2GO sequence table file as input.


usage: rename_goterms.py <infile.txt>
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	in_file_name = sys.argv[1]
	
goterms_file_name = "/Users/sinnafoch/Dropbox/Philipp/Genomes/000_important_stuff/goterms.txt" #laptop path!
goterms_file = open(goterms_file_name, "r")

in_file = open(in_file_name, "r")

goterms = {}
for line in goterms_file:
	goterms[line.strip().split("\t")[-1]] = line.strip().split("\t")[-2]


errors = []
i = 0
for line in in_file:
	try:
		i += 1
		name = line.split("\t")[0]
		terms = line.split("\t")[7]
		if ":" not in terms: #skip other lines
			continue
		terms = [term.strip(" ") for term in terms.split(";")]
		#print terms
		#time.sleep(1)
		terms = [item[2:] for item in terms]
		#print terms
		term_string = ""
		for term in terms:
			term_string += goterms[term] + ", "
		print name + "\t" + term_string[:-2]
		#not needed:
		#gene_id = "%06d" %(i,)
		#print gene_id + "\t" + term_string[:-2]
	except KeyError:
		sys.stderr.write("Problem with " + term +" in "+ name + "\n")
		errors.append(term)
		
sys.stderr.write("\nProblematic terms: " + str(set(errors)) + "\n")
sys.stderr.write("Total number: " + str(len(set(errors)))+"\n")