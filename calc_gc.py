#!/usr/bin/env python
# Philipp Resl, 100216
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_b

if len(sys.argv) < 2:
	print "calc_gc.py <infile>"
	quit()
else:
 	Filename = sys.argv[1]

window = 100
minlen = 500
stepsize = 0

File = open(Filename,"r")

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

listochunks = []
listogc = []
listonames = []
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
	#i += 1
	#if (i == 5): 
	#	break

maxlength = max(max(listochunks))
#print longest_seq
pdf = pdf_b.PdfPages("output.pdf")	
for name, chunk, gc in zip(listonames, listochunks, listogc):
	print "Plotting ", name
	with plt.style.context('fivethirtyeight'):
		fig = plt.figure()
		pl = fig.add_subplot(111)
		#pl.plot(chunk, gc, linestyle="None", marker=".", markerfacecolor="red",markersize=8)
		pl.scatter(chunk, gc, c=gc, s=300, alpha=0.6)
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
File.close()

	