#!/usr/bin/env python
# written by Philipp Resl
# 07.07.16

import sys
import re

info = """
Parses BLAST XML files (outfmt -5) created with piped output from parallel. 

Usage:
parse_blastxml.py <infile>

"""

if len(sys.argv) < 2:
	sys.stderr.write(info)
	quit()
else:
	file_name = sys.argv[1]

file = open(file_name, "r")
i = 1
excludes = ["<?xml version=\"1.0\"?>", "<!DOCTYPE","<BlastOutput>","<BlastOutput_program>","<BlastOutput_version>","<BlastOutput_reference>","<BlastOutput_db","<BlastOutput_query","<BlastOutput_param>","<Parameters>","<Parameters_","</Parameters>","</BlastOutput_param>","<BlastOutput_iterations>", "<BlastOutput_iterations>","</BlastOutput_iterations>","</BlastOutput>"]

def check_excludes(string):
	for exclude in excludes:
		if exclude in string:
			return 1
	return 0

head = file.readlines()[0:20]
for line in head:
	print line.strip("\n")
file.seek(0,0)
	
for line in file:
	if line == "\n": #get rid of empty lines
		continue
	line = line.strip("\n") #remove line feed
	if check_excludes(line) == 1: #check for lines to exclude
		continue	
	if "<Iteration_iter-num>" in line: #change iter number
		line = re.sub(r"\d+",str(i),line)
		i += 1
	print line

print "</BlastOutput_iterations>"
print "</BlastOutput>"
print

file.close()
		