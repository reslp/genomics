#!/usr/bin/env python
# Philipp Resl
# GENOMICS
# 09.07.16
import sys
import os
import zipfile
import re

Info = """
Formats the BLAST Output from Blast2Go into a single XML file for further use.
Expects that BLAST Output is kept as individual Zip files in one directory.

usage: blast2go_to_tab.py <directory>
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	dirname= sys.argv[1]

def calcmiss(seq1, seq2):
	return sum(seq1!=seq2 for seq1,seq2 in zip(seq1,seq2))

cwd = os.getcwd()
cwd += "/"+dirname
names_list = os.listdir(cwd)
names_list = [item for item in names_list if ".zip" in item] #remove non zip files
for zip_name in names_list:
	zf = zipfile.ZipFile(cwd+"/"+zip_name, "r")
	for xml_file in zf.infolist():
		ifile = zf.open(xml_file)
		for line in ifile.readlines():
			if "<query-title>" in line:
				qseqid = re.search('%s(.*)%s' % ("e>", "</"), line).group(1)
				continue
			if "<id>" in line:
				sseqid = re.search('%s(.*)%s' % ("d>", "</"), line).group(1)
				continue
			if "<accession>" in line:
				sacc = re.search('%s(.*)%s' % ("n>", "</"), line).group(1)
				continue
			if "<title>" in line:
				sseqid = re.search('%s(.*)%s' % ("e>", "</"), line).group(1)
				continue
			if "<taxid>" in line:
				taxid = re.search('%s(.*)%s' % ("d>", "</"), line).group(1)
				continue
			if "<sciname>" in line:
				sciname = re.search('%s(.*)%s' % ("e>", "</"), line).group(1)
				continue
			if "<bit-score>" in line:
				bitscore = re.search('%s(.*)%s' % ("e>", "</"), line).group(1)
				continue
			if "<score>" in line:
				btop = re.search('%s(.*)%s' % ("e>", "</"), line).group(1)
				continue
			if "<evalue>" in line:
				evalue = re.search('%s(.*)%s' % ("e>", "</"), line).group(1)
				continue
			if "<identity>" in line:
				nident = re.search('%s(.*)%s' % ("y>", "</"), line).group(1)
				continue
			if "<positive>" in line:
				pident = re.search('%s(.*)%s' % ("e>", "</"), line).group(1)
				continue
			if "<query-from>" in line:
				qstart = re.search('%s(.*)%s' % ("m>", "</"), line).group(1)
				continue
			if "<query-to>" in line:
				qend = re.search('%s(.*)%s' % ("o>", "</"), line).group(1)
				continue
			if "<hit-from>" in line:
				sstart = re.search('%s(.*)%s' % ("m>", "</"), line).group(1)
				continue
			if "<hit-to>" in line:
				send = re.search('%s(.*)%s' % ("o>", "</"), line).group(1)
				continue
			if "<align-len>" in line:
				length = re.search('%s(.*)%s' % ("n>", "</"), line).group(1)
				continue
			if "<gaps>" in line:
				gapopen = re.search('%s(.*)%s' % ("s>", "</"), line).group(1)
				continue
			if "<qseq>" in line:
				qseq = re.search('%s(.*)%s' % ("q>", "</"), line).group(1)
				continue
			if "<hseq>" in line:
				hseq = re.search('%s(.*)%s' % ("q>", "</"), line).group(1)
				continue
			if "</Hit>" in line:
				mismatch = calcmiss(qseq,hseq)
				print qseqid + "\t" + sseqid + "\t" + pident +"\t" + length +"\t" +str(mismatch) + "\t" + gapopen + "\t" + qstart + "\t" +qend +"\t" +sstart +"\t" +send +"\t" + evalue +"\t" +bitscore
				outline = ""
				continue

#standard blast column oder: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore



