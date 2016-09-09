#!/usr/bin/env python
# written by Philipp Resl
# 05.09.16

import sys
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_b
import numpy as np

info = """Usage:

gc_plot.py controlfile.txt


Control file should look like this:
/full/path/to/ggf3_file.gff3
/full/path/to/assembly_file.fasta
/full/path/to/megan_file.txt
<window size of sliding window (integer)>
<number of contigs to analyze (int)>
/path/to/output.pdf
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	file_name = sys.argv[1]
	
#change this to be dynamic

print "Reading control file..."
file = open(file_name, "r")
gff_file = file.readline().strip()
assemb_file = file.readline().strip()
megan_file = file.readline().strip()
window = int(file.readline().strip())
number = int(file.readline().strip())
outfile_name = file.readline().strip()
file.close()
print "Parameters from control file:"
print "GFF3 file:", gff_file
print "FASTA file:", assemb_file
print "MEGAN file:", megan_file
print
print "Window size:", str(window)
print "Number of contigs to analyze:", str(number)

#globals:
#window = 100
minlen = 5*window
stepsize = 0
#number=60
listochunks = []
listogc = []
listonames = []

def slide(seq):
	nchunks = 0
	gclist = []
	while len(seq) >= window:
		chunk = seq[:window]
		at = chunk.count("A") + chunk.count("T")
		gc = chunk.count("G") + chunk.count("C")
		gc_content = (gc/float(window))*100
		gclist.append(gc_content)
		seq = seq[window-stepsize:]	
		nchunks += 1
	chunklist = range(nchunks)
	return (chunklist, gclist)
#slide

def read_sequence_information():
	print "Opening fasta file:",assemb_file
	File = open(assemb_file,"r")
	i = 0
	longest_seq = 0
	for Sequence in SeqIO.parse(File, "fasta"):
		if len(Sequence.seq) < minlen:
			print "Skipping ", Sequence.name
			continue
		print "Reading ", Sequence.name
		if (len(Sequence.seq) > longest_seq):
			longest_seq = len(Sequence.seq)
		chunkn, gc_from_chunk = slide(Sequence.seq)
		listochunks.append(chunkn)
		listogc.append(gc_from_chunk)
		listonames.append(Sequence.name)
		i += 1
		if (i == number): 
			break
#read_sequence_information

	
def read_gene_information():
	
	i = 0
	print "Opening GFF file:", gff_file
	file = open(gff_file, "r")
	print "Reading file..."
	df = pd.read_csv(file,skiprows=1, error_bad_lines=False, sep="\t", header=None)
	names_list = list(set(df[0].tolist()))
	df = df.set_index(0)
	file.close()
	id_matches = ["gene"] #check if this correct
	genes_dict = {}
	gene_names_dict = {}
	print "Extraction gene positional information..."
	for name in names_list:
		df_sliced = df.loc[name]
		try:
			positions_x_b = df_sliced[df_sliced[2].isin(id_matches)][3].tolist()
			positions_x_e = df_sliced[df_sliced[2].isin(id_matches)][4].tolist()
			gene_names_dict[name] = list(df_sliced[df_sliced[2].isin(id_matches)][8].tolist())
			#genes_dict[name] = (sorted(list(set(positions_x_b)),key=int),sorted(list(set(positions_x_e)),key=int))
			genes_dict[name] = (list(positions_x_b),list(positions_x_e))
			i += 1
		except AttributeError: #this is crappy...
			print "Only one hit. Skipping", name
			continue
	for key in gene_names_dict.keys():
		gene_names_dict[key]=[item.split(";")[0].replace("ID=","") for item in gene_names_dict[key]]
	print "Genes retrieved from GFF file for", str(i), "contigs"
	return genes_dict, gene_names_dict
#read_gene_information

def read_megan_information():
	f = open(megan_file, "r")
	col_dict = {}
	for line in f:
		if "Bacteria" in line.split("\t")[1]:
			col = "red"
		elif "Ascomycota" in line.split("\t")[1]:
			col = "green"
		elif "Basidiomycota" in line.split("\t")[1]:
			col = "blue"
		else:
			col = "black"
		col_dict[line.split("\t")[0].replace("-mRNA-1", "")]=col
	f.close()
	return col_dict
#read_megan_information

def next_color():
	for item in x:
		yield item	


#gather information
genes, gene_names = read_gene_information()
read_sequence_information()
col_dict = read_megan_information()


#sys.exit(0)

#Here starts the plotting stuff
maxlength = max(max(listochunks))
pdf = pdf_b.PdfPages(outfile_name)	
for name, chunk, gc in zip(listonames, listochunks, listogc):
	print "Plotting ", name
	fig = plt.figure()
	with plt.style.context('fivethirtyeight'):
		pl = fig.add_subplot(111)
		#pl.plot(chunk, gc, linestyle="None", marker=".", markerfacecolor="red",markersize=8)
		pl.scatter(chunk, gc, c=gc, s=300, alpha=0.6)
		try:
			gene_x_b = genes[name][0]
			gene_x_b = map(float, gene_x_b)
			gene_x_b = [x / window for x in gene_x_b]
			gene_x_e = genes[name][1]
			gene_x_e = map(float, gene_x_e)
			gene_x_e = [x / window for x in gene_x_e]
			colors = iter([col_dict[y] for y in gene_names[name]])
			gene_y = [20, 25, 15] * (len(gene_x_e)/3) #this probably does not always work!
			for x1,x2,y1,y2 in zip (gene_x_b, gene_x_e, gene_y, gene_y):
				pl.plot([x1, x2], [y1, y2], 'k-', lw=6, c=colors.next())
		except KeyError:
			print "Key Error: Maybe", name, "does not have genes" 
		size = fig.get_size_inches()
		dpi = fig.get_dpi()
		fig.set_size_inches(size[0]*20, size[1]*2)
		#set length of x axis
		pl.set_ylim([0,100])
		pl.set_xlim([0,maxlength])
		# rescale axes to kilobases
		ticks=pl.get_xticks().tolist()
		ticks = [((tick * window) / 1000) for tick in ticks]
		pl.set_xticklabels(ticks, fontsize=30)
		pl.set_yticklabels(pl.get_yticks().tolist(), fontsize= 30)
		#pl.rc("text",fontsize= 22)
		pl.set_title(name, fontsize=30)
		pdf.savefig(fig)
pdf.close()





