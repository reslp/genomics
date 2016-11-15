#!/usr/bin/env python
# GENOMICS
# Philipp Resl
# last changed 14.11.16

import sys
Info = """
Filters codeml output generated with codeml.py (summary.txt) for genes with pvalues < 0.05 and min 1 site under selection

usage: filter_codeml.py <infile.txt>
"""

def analyze_outfile(name):
	out_file = open("out/results.out."+name+"_H1","r")
	BEB = 0
	nsites = 0
	for line in out_file:
		if line.startswith("Bayes Empirical Bayes (BEB)"):
			BEB=1
		if "*" in line and BEB==1:
			nsites += 1
		if line.startswith("\n") and BEB == 1:
			break
	out_file.close()
	return nsites

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	in_file_name = sys.argv[1]


sum_file = open(in_file_name,"r")
for line in sum_file:
	#print line.split()
	try:
		if float(line.split()[4]) < 0.05:
			sites=analyze_outfile(line.split()[0])
			if sites > 0:
				print line.strip()+"\t"+str(sites)+"sites"
	except ValueError:
		continue	
