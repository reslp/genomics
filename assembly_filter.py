#!/usr/bin/env python
from Bio import SeqIO
import numpy as np
from scipy import stats

assembly_file = "/Users/sinnafoch/Dropbox/Philipp/Genomes/00_data/Xylographa_parallela/01_genome/Xyl_par_genomeV1.fasta"
blast_file = "/Users/sinnafoch/Dropbox/Philipp/Genomes/00_data/Xylographa_parallela/06_blast/xyl_par_makerrun1_blastp_results_no_header.txt"
coverage_file = "/Users/sinnafoch/Dropbox/Philipp/Genomes/00_data/Xylographa_parallela/01_genome/coverage_per_contig_xylographa_v1.txt"

assembly = open(assembly_file,"r")
blast_results = open(blast_file,"r")
coverage_results = open(coverage_file,"r")

def calc_gc(sequence):
	GC = sequence.count("G")+sequence.count("C")
	return GC / float(len(sequence))

def count_blast(contig_id):
	blast_results.seek(0)
	bacteria = 0
	eukaryote = 0
	for result in blast_results:
		if contig_id in result and "Bacteria" in result:
			bacteria += 1
		if contig_id in result and "Eukaryota" in result:
			eukaryote +=1
	if bacteria == 0:
		return 0.0
	else:
		return float(bacteria) / float(eukaryote+bacteria)
		
def get_coverage(contig_id):
	coverage_results.seek(0)
	for result in coverage_results:
		if contig_id in result:
			items = [splits for splits in result.split(" ") if splits is not ""]
			return float(items[-1].strip())
				

contigs = list(SeqIO.parse(assembly, "fasta"))

i=0
#print "All contigs:"
#print "Name\tGCContent\tPro.ofBactHits\tConverage"

gc_list = []
bact_list = []
cov_list = []
id_list= []

for contig in contigs:
	gc = calc_gc(contig.seq)
	bact = count_blast(contig.id)
	cov = get_coverage(contig.id)
	gc_list.append(gc)
	bact_list.append(bact)
	cov_list.append(cov)
	id_list.append(contig.id)
	#i+=1
	#if i == 101:
	#	break

## select bad contigs

mean_gc = np.mean(gc_list[1:100])
sigma_gc = np.std(gc_list[1:100])
conv_gc = stats.norm.interval(0.95, loc=mean_gc, scale=sigma_gc) #does only work if data is normally distributed

mean_cov = np.mean(cov_list[1:100])
sigma_cov = np.std(cov_list[1:100])
conv_cov = stats.norm.interval(0.95, loc=mean_cov, scale=sigma_cov) #does only work if data is normally distributed

gc_variance = 0.08

print "Name\tGCContent\tPro.ofBactHits\tCoverage\tAnnotation"
good_contigs = []
for (gc, bact, cov, id) in zip(gc_list, bact_list, cov_list, id_list):
	append = ""
	if (cov < conv_cov[0] or cov > conv_cov[1]):
		append += "coverage "
	if gc > conv_gc[1] or gc < conv_gc[0]:
		append += "gc-content "
	if (bact >= 0.66):
		append += "blast "		
	print id + "\t" + str(gc) + "\t" + str(bact) + "\t" + str(cov) + "\t" + append	
	if append == "" or append == "coverage " or append=="gc-content " or append=="blast ":
		good_contigs.append(id)
	

print "GC List", len(gc_list)
print "Bact List", len(bact_list)
print "Cov List", len(cov_list)
print "ID List",len(id_list)
print "GC INterval:", conv_gc
print "Is GC Content normaly distributed:", stats.normaltest(gc_list[1:100])
print "Mean Cov:", conv_cov
print "Is Coverage Content normaly distributed:", stats.normaltest(cov_list[1:100])
print
print "List of good contigs:"
for contig in good_contigs:
	print contig
