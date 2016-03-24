#!/usr/bin/env python
#last edit: 22.03.2016

import sys
import pandas as pd #needs >0.16
from Bio import SeqIO

Info = """Finds orthologous sets of genes from website based OrthoMCL runs.
usage: ortho_select.py <control_file.txt> [-single]
	-single flag will split output into single locus files for orthologous groups (this may create many files!)
the control file should look like this:
species1
/path/to/species1_orthmcloutputfile
/path/to/species1_fasta_file
species2
/path/to/species2_orthmcloutputfile
/path/to/species2_fasta_file
...and so on
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	file_name = sys.argv[1]
	args = sys.argv[2:len(sys.argv)]

single = False
if ("-single" in args):
	print "Splitting sequences to single locus files"
	single = True

## read and import control file
print "Reading control file..."
control_file =  open(file_name, "r")
control_list = []
for name, ortho_path, fasta_path in zip(control_file, control_file, control_file):
	control_list.append((name.strip(), ortho_path.strip(), fasta_path.strip()))
	#print name.strip() 
	#print ortho_path.strip()
	#print fasta_path.strip()

## import ortho_files and reduce to orthologous groups that are present in all species
ortho_groups_list = []
for species in control_list:
	print "\nExtracting Orthogroups for", species[0]
	ortho_file = pd.read_csv(species[1],sep="\t", header=None)
	total = list(ortho_file[1])
	print "Total number of records:", len(total)
	unique = list(set(ortho_file[1]))
	print "Number of unique groups:", len(unique)
	double = set([x for x in list(ortho_file[1]) if list(ortho_file[1]).count(x) > 1])
	print "Number of doubled groups:", len(double)
	double2 = [x for x in unique if x not in double]
	print "Number of unique groups (without doubles):", len(double2)
	#x = ortho_file[ortho_file[1] in double2	
	ortho_groups_list.append(double2)


## extract Sequence IDs
orthogroups = list(set(ortho_groups_list[0]).intersection(*ortho_groups_list))
print "\nNumber of shared Orthogroups between all species:", len(orthogroups)
id_list = []
orthos_list = []
for species in control_list:
	ortho_file = pd.read_csv(species[1],sep="\t", header=None)
	ids = ortho_file[ortho_file[1].isin(orthogroups)]
	id_list.append(list(ids[0]))
	orthos_list.append(list(ids[1]))

## extract sequences from fasta files and create new files for each locus -single flag
if single == True:
	print "\nWriting records to single locus files:"
	i=0
	j=1
	n = 500 #number of genes to extract, if you want to decrease runtime
	for orthogroup in orthogroups:
		print str(j), " Writing file for", orthogroup
		outfile = open(orthogroup +" _"+ str(len(control_list)) + "_species.fasta","w")
		for species in control_list:
			seqfile = open(species[2], "r")
			seqs_list = list(SeqIO.parse(seqfile, "fasta"))
			for sequence in seqs_list:
				index = orthos_list[i].index(orthogroup)
				if sequence.id == id_list[i][index]:
					sequence.id = orthos_list[i][index] + "_" +species[0] + "_" + sequence.id
					SeqIO.write(sequence, outfile, "fasta")	
			seqfile.close()
			i += 1
		i = 0
		j += 1
		outfile.close()
		if n == j:
			break
		

## creates files for each species
else:	
	print "\nWriting records to species files:"
	i = 0

	for species in control_list:
		print " Extracting sequences of shared Orthogroups for:", species[0]
		seqfile = open(species[2], "r")
		seqs_list = list(SeqIO.parse(seqfile, "fasta"))
		outfile = open(species[0] + "_orthologous_protein.fas", "w")
		for sequence in seqs_list:
			if sequence.id in id_list[i]:
				index = id_list[i].index(sequence.id)
				sequence.id = orthos_list[i][index] + "_" +species[0] + "_" + sequence.id
				SeqIO.write(sequence, outfile, "fasta")
		outfile.close()	
		i += 1
		seqfile.close()
