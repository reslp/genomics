#!/usr/bin/env python
# this script converts GFF3 with MAKER predictions to GTF for use with HISAT2

import sys

Info = """
This script converts GFF3 files as created by MKER2 to GTF files for use with HISAT2/Springtie

example: gff3_to_gtf.py <infile > <outfile>
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	gff_file_name = sys.argv[1]

gff_file = open(gff_file_name, "r")

for line in gff_file:
	if len(line.split("\t")) != 9:
		continue
	else:
		elements = line.split("\t")
		if elements[2] == "gene":
			gene_id = elements[8].split(";")[0].replace("ID=", "gene_id \"")+"\""
			print elements[0] + "\t" + elements[1]+ "\t" + elements[2]+ "\t" + elements[3]+ "\t" + elements[4]+ "\t" + elements[5]+ "\t" + elements[6]+ "\t" + elements[7]+ "\t" +gene_id + "; "
		elif elements[2] in ["CDS","exon"]:
			transcript_id = elements[8].split(";")[1].replace("Parent=", "transcript_id \"").strip()+"\""
			print elements[0] + "\t" + elements[1]+ "\t" + elements[2]+ "\t" + elements[3]+ "\t" + elements[4]+ "\t" + elements[5]+ "\t" + elements[6]+ "\t" + elements[7]+ "\t" +gene_id + "; " + transcript_id +";"
		else:
			continue
			
				