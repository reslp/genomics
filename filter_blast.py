#!/usr/bin/env python
import sys
import pandas as pd #needs >0.16

Info = """
Filters multiple blast results for best e-value. sorts by query name (default=qseqid) and evalue (default=evalue).
BLAST results need to be in tabular format and with column names. see default values above.


usage: filter_blast.py <blast_results.txt>
"""

if len(sys.argv) < 2:
	sys.stderr.write(Info)
	quit()
else:
	file_name = sys.argv[1]

sequence_id="qseqid"
evalue_id="evalue"
	
Database = pd.read_csv(file_name,sep="\t")
names = list(set(Database[sequence_id].tolist())) #remove duplicates from list of names
new_db= []
for name in names:
	db_sub = Database.loc[Database[sequence_id] == name]
	evalue = db_sub[evalue_id].min()
	small_db = pd.DataFrame(db_sub.loc[db_sub[evalue_id] == evalue])
	if (len(small_db) > 1): #if multiple hits have same e-value
		small_db = small_db.sample(n=1) # get random row when e-values are the same
	if small_db["sskingdoms"].to_string(index=False) == "Eukaryota":
		new_db.append(small_db)
		
pd.concat(new_db).to_csv(sys.stdout,sep="\t", index=False)

	
