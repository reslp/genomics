#!/usr/bin/env python
import sys
import pandas as pd #needs >0.16

Info = """
Filters ORTHOMCL output for specific Orthogroups. 


usage: filter_orthomcl.py <orthmcl_results.txt> <query_file.txt>
"""

if len(sys.argv) < 3:
	sys.stderr.write(Info)
	quit()
else:
	file_name = sys.argv[1]
	query = sys.argv[2]
	
Database = pd.read_csv(file_name,sep="\t")
#names = list(set(Database[sequence_id].tolist())) #remove duplicates from list of names
new_db= []

file = open(query,"r")

for query in file:
	query = query.strip()
	for iter, item in Database.iterrows():
		if query == item[1]:
			print item[0]


file.close()
#pd.concat(new_db).to_csv(sys.stdout,sep="\t", index=False)

	
