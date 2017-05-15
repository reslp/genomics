#!/usr/bin/env python
# Philipp Resl
# GENOMICS
# 11.01.17
import sys
from collections import Counter

Info = """
THis script summarizes the occurence of GO-Terms. For each GO Term in a list of terms it gives the occurence 
frequency in the complete list.

usage: summarize_go.py <infile>


Example for Format of Infile:

C:cellular_component 
F:molecular_function 
P:biological_process 
F:structural constituent of ribosome 
P:translation 
C:mitochondrion 
C:ribosome 
C:cytosolic small ribosomal subunit 
F:proton-transporting ATPase activity, rotational mechanism 
C:nuclear pore 
F:structural constituent of nuclear pore 
P:nucleocytoplasmic transport 
P:aerobic respiration 
P:ATP synthesis coupled proton transport 
P:hydrogen ion transmembrane transport 
C:cytosolic large ribosomal subunit 

"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	filename= sys.argv[1]
	

file = open(filename, "r")

goterms = []

for line in file:
	goterms.append(line.strip())
	
godict = Counter(goterms)

	
for item in sorted(godict, key=godict.get, reverse=True): #godict.keys():
	print item, ": ", godict[item]

line = []
for item in goterms:
	item = item.replace(":", " ")
	line.append(item.split())

flatten = [item for sublist in line for item in sublist]
print flatten

worddict = {}
for item in flatten:
	#print item + " " + str(flatten.count(item))
	worddict[item] = flatten.count(item)

for item in sorted(worddict, key=worddict.get, reverse=False): #godict.keys():
	print item, ": ", worddict[item]
