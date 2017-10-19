#!/usr/bin/env python
# written by Philipp Resl

import numpy as np
from scipy import misc
from PIL import Image
import sys


Info = """
Creates a visual image from a genome assembly in PNG format.

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
for line in file:
	if line.startswith(">")==False:
		sequence += line.strip()

print "Assembly length: %s" % str(len(sequence))
x = len(sequence) / 10000
y = 10001

im = Image.new("RGB",(x,y))
data = np.zeros((y,x,3), dtype=np.uint8)

xx = 0
yy = 0
for letter in sequence:
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
	if xx == x: #to create image at defined width
		yy +=1
		xx = 0
	if yy == y:
		break

img = misc.toimage(data) 
img.save(out_file_name)
