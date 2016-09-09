#!/usr/bin/env python
#last edit: 09.09.2016

import sys
import pandas as pd #needs >0.16
from Bio import SeqIO
from operator import itemgetter # for sorting list of tuples

Info = """Finds orthologous sets of genes from website based OrthoMCL runs.
usage: ortho_select.py <control_file.txt> [-single] [-invert] -ava
	-single flag will split output into single locus files for orthologous groups (this may create many files!)
	-invert flag will only output sequences that are not present in other species. -single will not work with -invert
	-ava will read controlfile for ortholog data created by all vs. all blast
the control file should look like this:
species1
/path/to/species1_orthmcloutputfile
/path/to/species1_fasta_file
species2
/path/to/species2_orthmcloutputfile
/path/to/species2_fasta_file
...and so on

control file for -ava mode should look like this:
/path/to/orthofile # orthofile needs to have column names identical to species names in control file
species1
/path/to/species1_fasta_file
species2
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
invert = False
ava = False
if "-single" in args:
	print "Splitting sequences to single locus files"
	single = True
if "-invert" in args:
	print "Inverted selection. Will output only unique orthogroups"
	invert = True
if "-ava" in args:
	print "Mode is set to interpret custom all vs. all blast files created with get_best_hit.py and filter_reciprocal.py"
	ava = True

## read and import control file
print "Reading control file..."
if ava == True:
	if invert == False:
		control_file = open(file_name, "r")
		ortho_file_name = control_file.readline().strip()
		control_list = []
		for name, fasta_path in zip(control_file, control_file):
			control_list.append((name.strip(), fasta_path.strip()))
		control_file.close()
		print "File with ortholog information: ", ortho_file_name
		ortho_file = pd.read_csv(ortho_file_name, sep ="\t", header=None, error_bad_lines=True)
		ortho_file.columns = ortho_file.iloc[0] 
		
		#print control_list
		print "Total number of records: ", len(ortho_file)-1
		if single == True:
			seq_dict = {}
			for species in control_list:
				genelist = list(ortho_file[species[0]])
				#print len(genelist)
				del genelist[0]
				infile = open(species[1],"r")
				seqs_list = list(SeqIO.parse(infile, "fasta"))
				#outfile = open(species[0]+"_orthologs.fas","w")
				ortholist= []
				for gene in genelist:
					j = 0
					for sequence in seqs_list:
						if gene in species[0] + "_" + sequence.id: #this is a workaround because SeqRecords are edited in place!
							#print gene, sequence.id
							#print sequence.id
							id = ">"+species[0] + "_" + sequence.id+"\n"
							seq = str(sequence.seq) +"\n"
							ortholist.append(id+seq)
							j = 1
					if j == 0:
						print "Problem:", gene, "not found"
				print str(len(ortholist)), "genes for species", species[0], "extracted"
				seq_dict[species] = ortholist
				infile.close()

			print "Writing records to single locus files:"
			for i in range(len(ortho_file)-1):
				file = open(str(i+1)+"_ortholog_"+str(len(control_list))+"species.fas", "w")
				print "Writing " + str(i+1)+"_ortholog_"+str(len(control_list))+"species.fas"
				keys = seq_dict.keys()
				keys = sorted(keys,key=itemgetter(0))
				for key in keys:
					file.write(seq_dict[key][i])
				file.close()			
		else: #single = False
			for species in control_list:
				print "Extracting sequences for ", species[0]
				genelist = list(ortho_file[species[0]])
				del genelist[0]
				infile = open(species[1],"r")
				seqs_list = list(SeqIO.parse(infile, "fasta"))
				outfile = open(species[0]+"_orthologs.fas","w")
				for gene in genelist:
					for sequence in seqs_list:
						if gene in species[0] + "_" + sequence.id: #this is a workaround because SeqRecords are edited in place!
							id = ">"+species[0] + "_" + sequence.id+"\n"
							seq = str(sequence.seq) +"\n"
							outfile.write(id)
							outfile.write(seq)	
				outfile.close()
				infile.close()
	else: #invert = True
		print "Use other script! to extract unique genes from ava."
		
			
	
else:
	control_file =  open(file_name, "r")
	control_list = []
	for name, ortho_path, fasta_path in zip(control_file, control_file, control_file):
		control_list.append((name.strip(), ortho_path.strip(), fasta_path.strip()))
	control_file.close()


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
	if (invert == False):
		orthogroups = list(set(ortho_groups_list[0]).intersection(*ortho_groups_list))
		print "\nNumber of shared Orthogroups between all species:", len(orthogroups)
		id_list = []
		orthos_list = []
		for species in control_list:
			ortho_file = pd.read_csv(species[1],sep="\t", header=None)
			ids = ortho_file[ortho_file[1].isin(orthogroups)]
			id_list.append(list(ids[0]))
			orthos_list.append(list(ids[1]))
	if (invert == True):
		orthos_list = []
		id_list = []
		#print len(ortho_groups_list[0])
		#print len(ortho_groups_list)
		temp_list_other_species = ortho_groups_list
		for i in range(0,len(ortho_groups_list)):
			ortho_file = pd.read_csv(control_list[i][1],sep="\t", header=None)
			temp_list_species = ortho_groups_list[i]
			temp_list_other_species = list(ortho_groups_list)
			temp_list_other_species.pop(i)
			temp_list_other_species = [item for sublist in temp_list_other_species for item in sublist]
			orthos_list.append([item for item in temp_list_species if item not in temp_list_other_species])
			ids = ortho_file[ortho_file[1].isin(orthos_list[i])]
			id_list.append(list(ids[0]))
			print "\nNumber of unique Orthogroups for species: ", control_list[i][0], len(orthos_list[i]), len(ids[0])

	## extract sequences from fasta files and create new files for each locus -single flag
	if single == True and invert == False:
		print "\nWriting records to single locus files:"
		i=0
		j=1
		n = 10000 # number of genes to extract, if you want to decrease runtime
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
			if invert == False:
				print "Extracting sequences of shared Orthogroups for:", species[0]
				seqfile = open(species[2], "r")
				seqs_list = list(SeqIO.parse(seqfile, "fasta"))
				outfile = open(species[0] + "_orthologous_protein.fas", "w")
				for sequence in seqs_list:
					if sequence.id in id_list[i]:
						index = id_list[i].index(sequence.id)
						sequence.id = orthos_list[i][index] + "_" +species[0] + "_" + sequence.id
						SeqIO.write(sequence, outfile, "fasta")
			else:
				print "Extracting sequences of unique Orthogroups for:", species[0]
				seqfile = open(species[2], "r")
				seqs_list = list(SeqIO.parse(seqfile, "fasta"))
				outfile = open(species[0] + "_unique_orthologous_protein.fas", "w")
				for sequence in seqs_list:
					if sequence.id in id_list[i]:
						index = id_list[i].index(sequence.id)
						#sequence.id = orthos_list[i][index] + "_" +species[0] + "_" + sequence.id
						SeqIO.write(sequence, outfile, "fasta")
			i += 1
	seqfile.close()
	outfile.close()	