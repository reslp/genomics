#!/usr/bin/env python
#last edit: 22.03.2016
#biopython needs to be 1.66!
from Bio import SeqIO
from Bio import AlignIO
from Bio import codonalign
from Bio.Alphabet import IUPAC
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.Align.Applications import MafftCommandline
from StringIO import StringIO
import numpy as np
from operator import truediv
from os import listdir
import sys

## here import all the single locus protein files
def listdir_nohidden(path):
    for f in listdir(path):
        if not f.startswith('.'):
            yield f
            
protein_file_names = listdir_nohidden("proteins")
dna_file_names = listdir_nohidden("transcripts")
#print protein_file_names
#print dna_file_names

## function to create alignment for individual files
def align_next_file (name):
	mafft_cmd = MafftCommandline("/usr/local/bin/mafft", input="proteins/"+name)
	stdout, stderr = mafft_cmd()
	alignment = AlignIO.read(StringIO(stdout), "fasta", alphabet=IUPAC.protein)
	return alignment
	
## function to load individual dna files
def import_next_file (name):
	fn = "transcripts/"+name
	file= open(fn,"r")
	nucl = SeqIO.parse(file,"fasta", alphabet=IUPAC.IUPACUnambiguousDNA())
	return nucl

#print "Gene\tdN \tdS\tdN/DS"
i = 0
n = 1000
print "gene\txylogr\ttrape\tagyri\tgraphi\tlambie\tcladon\txantho"
outfile = open("output.txt", "a")
outfile.write("gene\txylogr\ttrape\tagyri\tgraphi\tlambie\tcladon\txantho\n")
outfile.close()

for prot_file, dna_file in zip(protein_file_names, dna_file_names):
	outfile = open("output.txt", "a")
	prot = align_next_file(prot_file)
	dna = import_next_file(dna_file)
	codon = codonalign.build(prot, dna, alphabet=codonalign.default_codon_alphabet)
	#dN, dS = cal_dn_ds(codon[0], codon[1], method='NG86')
	dn_matrix, ds_matrix = codon.get_dn_ds_matrix(method="ML")
	length = len(sorted(dn_matrix.matrix,key=len,reverse=True)[0]) #although this is always the same depending on the number of species
	dn_mat = np.array([xi+[0.0]*(length-len(xi)) for xi in dn_matrix.matrix])
	ds_mat = np.array([xi+[0.0]*(length-len(xi)) for xi in ds_matrix.matrix])
	sum_dn_col = list(np.sum(dn_mat, axis=0))
	sum_dn_row = list(np.sum(dn_mat, axis=1))
	sum_dn_col = [x/length for x in sum_dn_col]
	sum_dn_row = [x/length for x in sum_dn_row]
	sum_ds_col = list(np.sum(ds_mat, axis=0))
	sum_ds_row = list(np.sum(ds_mat, axis=1))
	sum_ds_col = [x/length for x in sum_ds_col]
	sum_ds_row = [x/length for x in sum_ds_row]
	marg_dn = [x+y for x,y in zip(sum_dn_col,sum_dn_row)]
	marg_ds = [x+y for x,y in zip(sum_ds_col,sum_ds_row)]
	dnds = [x/y for x,y in zip(marg_dn,marg_ds)]
	output = prot_file[0:10] +"\t"
	for value in dnds:
		output += str("%.4f" % value) + "\t"
	print output
	outfile.write(output+"\n")
	#print prot_file,"\t", dN,"\t",dS,"\t", dN/dS
	outfile.close()
	i += 1
	if i == n:
		break
sys.exit(0)


#codon = codonalign.build(alignment, nucl, alphabet=codonalign.default_codon_alphabet)
#dN, dS = cal_dn_ds(codon[0], codon[1], method='NG86')
#print dN, dS
#print dN/dS
#dn_matrix, ds_matrix = codon.get_dn_ds_matrix(method="NG86")
#print dir(dn_matrix)
#print dn_matrix.matrix /ds_matrix.matrix
#a = [item for sublist in ds_matrix.matrix for item in sublist]
#b =[item for sublist in dn_matrix.matrix for item in sublist]
#c = map(truediv, dn_matrix.matrix,ds_matrix.matrix)
#print a/b

#print vars(dn_matrix)
