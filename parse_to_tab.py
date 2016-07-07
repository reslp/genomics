#!/usr/bin/env python
# written by Philipp Resl
# 07.07.16

import sys
from Bio import SearchIO

info = """
Creates tabular format from XML blast output file.

Usage:
parse_to_tab.py <infile>

"""

if len(sys.argv) < 2:
	sys.stderr.write(info)
	quit()
else:
	file_name = sys.argv[1]
	
file = open(file_name,"r")

blast_results = SearchIO.parse(file, 'blast-xml')
for result in blast_results:
	print result
	#SearchIO.write(result, sys.stdout, 'blast-tab')