#!/usr/bin/env python

import argparse
#import glob
import sys
import os
from Bio import SeqIO
import pandas as pd
import numpy as np

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="extract_annotation_for_genes.py", description = """This script extracts annotation data from funannotate annotation result files for a set of genes""", epilog = """written by Philipp Resl""")
pars.add_argument('-i', dest="i", required=True, help="Path to the funnannotate annotation file")
pars.add_argument('-g', dest="g", required=True, help="Path to FASTA file with sequences")
pars.add_argument('-t', dest="t", required=True, help="Type of annotation to be extracted. Can be: Product, BUSCO, PFAM, InterPro, EggNog, COG, GO Terms, Secreted, Membrane, Protease, CAZyme")


args=pars.parse_args()

valid_annotations = ["Name", "Product", "BUSCO", "PFAM", "InterPro", "EggNog", "COG", "GO Terms", "Secreted", "Membrane", "Protease", "CAZyme"]
if args.t not in valid_annotations:
	print("Wrong argument -t", args.t ,"supplied. Check again.", file=sys.stderr)
	sys.exit(1)

data = pd.read_csv(args.i, sep="\t", header=0)
data = data.set_index("GeneID")

genes = list(SeqIO.parse(args.g, "fasta"))
idlist = [gene.id.split("-")[0] for gene in genes]

#print(data)

for id in idlist:
	#print(id, args.t)
	annotations = str(data.loc[id, args.t])
	if annotations == "nan": #remove missing data
		annotations = ""
	print(id, annotations)
