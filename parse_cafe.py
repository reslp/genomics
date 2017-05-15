#!/usr/bin/env python
# written by Philipp Resl
# 01.27.17

import sys
import re
from collections import defaultdict

info = """
Parses output from cafe 3.0 to make it more readable

Usage:
parse_cafe.py <infile> [-ortho]

"""

if len(sys.argv) < 2:
	sys.stderr.write(info)
	quit()
else:
	file_name = sys.argv[1]
	flags = sys.argv[2]

if "-ortho" in flags:
	cafe_log_file = open(file_name, "r")
	key_file = open("/Users/sinnafoch/Dropbox/Philipp/Genomes/14_orthogroup_expansion/all_cafe.txt", "r")
	key_dict = {}
	for line in key_file:
		key = line.strip().split("\t")[1]
		value = ""
		key_dict[key] = value
else:
	cafe_log_file = open(file_name, "r")
	key_file = open("/Users/sinnafoch/Dropbox/Philipp/Genomes/000_important_stuff/pfam_all_description.txt", "r")
	key_dict = {}
	for line in key_file:
		#print line
		key, value = line.strip().split("\t")
		key_dict[key] = value

#print key_dict	
results_dict = {}
for line in cafe_log_file:
	if "(node ID, node ID):" in line:
		node_ids = line.split("(node ID, node ID): ")[1].strip()
		node_list = filter(None,re.split("[,() ]", node_ids))
		#print len(node_list)
	if line.startswith("PF"):
		line = line.strip()
		family_id = line.split("\t")[0]
		tree = line.split("\t")[1]
		global_pvalue = line.split("\t")[2]
		viterbi_pvalues = line.split("\t")[3]
		results_dict[family_id] = (tree, global_pvalue, viterbi_pvalues)
	if line.startswith("OG"):	
		line = line.strip()
		family_id = line.split("\t")[0]
		tree = line.split("\t")[1]
		global_pvalue = line.split("\t")[2]
		viterbi_pvalues = line.split("\t")[3]
		results_dict[family_id] = (tree, global_pvalue, viterbi_pvalues)

sign_families_at_node = defaultdict(list)
total_sign_count = 0
for family in results_dict.keys():
	if float(results_dict[family][1]) < 0.01:
		total_sign_count += 1
		#print family, results_dict[family][2]
		vit_p_list = [item.strip(")").strip("(") for item in results_dict[family][2].split(",")]
		#print len(vit_p_list)
		for i in range(0,len(vit_p_list)):
			#print i
			#print vit_p_list[i]
			if float(vit_p_list[i]) < 0.01:
				#print "check"
				sign_families_at_node[node_list[i]].append(family)

for node in sign_families_at_node.keys():
	print "Node", node
	print "Total number:", len(sign_families_at_node[node])
	for item in sign_families_at_node[node]:
		print item +"\t"+ key_dict[item]
	#print sign_families_at_node[node]
		
print "Total significant families:", total_sign_count
