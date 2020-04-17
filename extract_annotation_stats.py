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

pars = argparse.ArgumentParser(prog="extract_annotation_stats.py", description = """This script extracts annotation data from funannotate annotation result files.""", epilog = """written by Philipp Resl""")
pars.add_argument('-i', dest="i", required=True, help="Path to the funnannotate annotation file")
args=pars.parse_args()


data = pd.read_csv(args.i, sep="\t", header=0)

#print(data)

print("No. of predicted genes:\t", len(data["GeneID"]))
no_of_nans = data[["Secreted", "PFAM", "InterPro", "EggNog", "GO Terms", "CAZyme", "Membrane", "Protease"]].isnull().sum(axis=1).tolist()
print("No. of Genes without annotation:\t", no_of_nans.count(8))
sec_annotations = data["Secreted"].dropna().tolist()
print("No. of secreted genes:\t", len(sec_annotations))


pfam_annotations = data["PFAM"].dropna().tolist()
print("No. of genes with 1+ PFAM annotations:\t", len(pfam_annotations))
pfams = [ann for anns in pfam_annotations for ann in anns.split(";")]
print("Total no. of PFAM annotations:\t", len(pfams))
print("Unique no. of PFAM domains:\t", len(set(pfams)))

inter_annotations = data["InterPro"].dropna().tolist()
print("No. of genes with 1+ InterPro annotations:\t", len(inter_annotations))
interp = [ann for anns in inter_annotations for ann in anns.split(";")]
print("Total no. of InterPro annotations:\t", len(interp))
print("Unique no. of InterPro domains:\t", len(set(interp)))

egg_annotations = data["EggNog"].dropna().tolist()
print("No. of genes with EggNog Orthologs:\t", len(egg_annotations))
eggn = [ann for anns in egg_annotations for ann in anns.split(";")]
print("Total no. of EggNog Orthologs:\t", len(eggn))
print("Unique no. of EggNog Orthologs:\t", len(set(eggn)))

go_annotations = data["GO Terms"].dropna().tolist()
print("No. of genes with 1+ GO Term:\t", len(go_annotations))
gon = [ann for anns in go_annotations for ann in anns.split(";")]
print("Total no. of GO Terms:\t", len(gon))
print("Unique no. of GO Terms:\t", len(set(gon)))

cazy_annotations = data["CAZyme"].dropna().tolist()
#print(cazy_annotations)
print("No. of genes CAZymes:\t", len(cazy_annotations))
cazyn = [ann for anns in cazy_annotations for ann in anns.split(";")]
print("Total no. of CAZymes:\t", len(cazyn))
print("Unique no. of CAZymes:\t", len(set(cazyn)))

mem_annotations = data["Membrane"].dropna().tolist()
#print(mem_annotations)
print("No. of Membrane Proteins:\t", len(mem_annotations))
memn = [ann for anns in mem_annotations for ann in anns.split(";")]
print("Total no. of Membrane Proteins:\t", len(memn))
print("Unique no. of Membrane Proteins:\t", len(set(memn)))

prot_annotations = data["Protease"].dropna().tolist()
#print(prot_annotations)
print("No. of Proteases:\t", len(prot_annotations))
protn = [ann for anns in prot_annotations for ann in anns.split(";")]
print("Total no. of Proteases:\t", len(protn))
print("Unique no. of Proteases:\t", len(set(protn)))




#data_subs = data[["GeneID", "Contig", "Start", "Stop"]]
#data_subs.to_csv("genes_map.csv", sep="\t", header=False, index=False)