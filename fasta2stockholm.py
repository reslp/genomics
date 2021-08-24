#!/bin/python3
# written by Philipp Resl
import sys
from Bio import AlignIO

name = sys.argv[1]
#align = AlignIO.read(name, "fasta")
AlignIO.convert(name, "fasta", sys.stdout, "stockholm", alphabet=None)