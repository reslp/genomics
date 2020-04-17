#!/usr/bin/env python

import argparse
#import glob
import sys
import os
#from Bio import SeqIO
import pandas as pd
import numpy as np

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="create_gene_universe_file.py", description = """This script extracts GO annotation data from a funannotate pseudo GFF3 file""", epilog = """written by Philipp Resl""")
pars.add_argument('-i', dest="i", required=True, help="Path to the funnannotate annotation pseudo GFF file")
args=pars.parse_args()


data = pd.read_csv(args.i, sep="\t", header=0)

#go_annotations = data["GO Terms"].tolist()

geneids = data["GeneID"].tolist()

for gene in geneids:
	go_list = data.loc[data["GeneID"] == gene, "GO Terms"].fillna("").tolist()
	go_list = "".join(go_list)
	go_list = go_list.replace(";", ",")
	print(gene,"\t", go_list)
	
#print("No. of genes with 1+ GO Term:\t", len(go_annotations))
#gon = [ann for anns in go_annotations for ann in anns.split(";")]

#data
#print("Total no. of GO Terms:\t", len(gon))
#print("Unique no. of GO Terms:\t", len(set(gon)))






#data_subs = data[["GeneID", "Contig", "Start", "Stop"]]
#data_subs.to_csv("genes_map.csv", sep="\t", header=False, index=False)