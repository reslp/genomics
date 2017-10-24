#!/usr/bin/env python
# written by Philipp Resl

import numpy as np
from scipy import misc
import sys


Info = """
Creates a visual image from a genome assembly in PNG format.
Due to constraints in the image creating process assemblies >100MB will be split up in multiple files.

usage: $ genome_vix.py assembly.fasta outfile.png
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	assembly_file_name = sys.argv[1]
	out_file_name = sys.argv[2]

file = open(assembly_file_name, "r")

sequence = ""
chunk_size=100000000 #100 Megabases

for line in file:
	if line.startswith(">")==False:
		sequence += line.strip()

print "Assembly length: %s" % str(len(sequence))

if len(sequence) > chunk_size:
	print "The assembly is large (>100MB). Will create multiple files"

# the size of the image
x = 10000
y = 10000

data = np.zeros((y,x,3), dtype=np.uint8)

xx = 0
yy = 0
i = 0
img_count = 0
for letter in sequence:
	i += 1
	if letter =="A":
		data[yy,xx] = [254,0,0] #red
	elif letter == "G":
		data[yy,xx] = [0,0,254] #blue
	elif letter == "C":
		data[yy,xx] = [0,254,0] #green
	elif letter == "T":
		data[yy,xx] = [254,254,0] #yellow
	else:
		data[yy,xx] = [128,128,128] #grey
	xx += 1
	if xx == x: #to create image with defined width
		yy +=1
		xx = 0
	if yy == y or i == len(sequence):
		img_count += 1
		img = misc.toimage(data[0:yy,0:x])
		out_file = str(img_count) + "_" + out_file_name 
		img.save(out_file)
		xx = 0
		yy = 0
	


