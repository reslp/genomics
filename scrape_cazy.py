#!/usr/bin/env python

import argparse
import requests
import urllib.request
import time
from bs4 import BeautifulSoup
import sys

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="scrape_cazy.py", description = """This script downloads and parses CAZyme information from cazy.org""", epilog = """written by Philipp Resl""")
pars.add_argument('-f', dest="fam", required=True, help="List (delimited by commas) of CAZyme families")
args=pars.parse_args()
cazy_list = args.fam.split(",")


def get_cazy_table(cazy):
	cazy_data = ""
	url = "http://www.cazy.org/%s_characterized.html" % cazy
	response = requests.get(url)
	soup = BeautifulSoup(response.text, "html.parser")
	# First, get the number of charcacterized Enzymes, this is important for deciding if there are multiple pages.
	# This happens if the number if > 100
	empty = ""
	if soup.find(id="line_actif") == None:
		print("Something wrong with", url, file=sys.stderr)
		return ""
	if soup.find(id="line_actif").get_text().split("(")[-1].split(")"[0]) != empty:
		number = int(soup.find(id="line_actif").get_text().split("(")[-1].split(")")[0])
		print("Found a total number of", number, "cazymes in family", cazy, file=sys.stderr)

	#print(soup.find(id="line_actif").get_text())
	#print(soup.prettify())

	for page in range(0, number, 100):
		if page == 0:
			tables = soup.findChildren("table")
			cazy_table = tables[1]
			rows = cazy_table.findChildren("tr")
			rows = rows[1:]
			prefix = cazy + "\t"
			title = True
			for row in rows:
				#print(row)
				if "line_titre" in str(row) and title == True:
					title = False
				elif "line_titre" in str(row) and title == False:
					continue
				if row.find("td", {"class": "separateur1"}) != None:
					if "Top" in row.find("td", {"class": "separateur1"}).get_text():
						continue
					else:
						prefix += row.find("td", {"class": "separateur1"}).get_text()
						prefix += "\t"
						continue
					#print(row.find("td", {"class": "separateur1"}).get_text())
				columns = row.findChildren("td")
				string = prefix
				for column in columns:
					if "Top" in column.get_text():
						continue
					else:
						string += column.get_text(" ").strip()
						string += "\t"
				string += "\n"
				cazy_data += string
		else:
			add_to_url = "?debut_FUNC=%s#pagination_FUNC" % str(page)
			response = requests.get(url + add_to_url)
			soup = BeautifulSoup(response.text, "html.parser")
			tables = soup.findChildren("table")

			cazy_table = tables[1]
			rows = cazy_table.findChildren("tr")
			rows = rows[1:] # skip first two lines which contain the headers
			prefix = cazy + "\t"
			for row in rows:
				if "line_titre" in str(row):
					continue
				if row.find("td", {"class": "separateur1"}) != None:
					if "Top" in row.find("td", {"class": "separateur1"}).get_text():
						continue
					else:
						prefix += row.find("td", {"class": "separateur1"}).get_text()
						prefix += "\t"
						continue
					#print(row.find("td", {"class": "separateur1"}).get_text())
				columns = row.findChildren("td")
				string = prefix
				for column in columns:
					if "Top" in column.get_text():
						continue
					else:
						string += column.get_text(" ").strip()
						string += "\t"
				string += "\n"
				cazy_data += string
	return cazy_data

for arg in cazy_list:
	cazy = arg
	filename = cazy + "_characterized.txt"
	outfile = open(filename, "w")
	print(get_cazy_table(cazy), file = outfile)
	outfile.close()
